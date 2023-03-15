from collections import OrderedDict

import numpy as np
import pytorch_lightning as pl
import torch
import torch.nn.functional as F
from torch import nn


# ------------------------------------------------------------------------------
# ConvNext1D Components
# ------------------------------------------------------------------------------
class GlobalResponseNormalization(nn.Module):
    """
    Adapted from https://github.com/facebookresearch/ConvNeXt-V2

    Inputs are expected to be of shape (batch_size, channels, spatial_dim),
    just as with conv1d layers.
    """

    def __init__(self, dim):
        super(GlobalResponseNormalization, self).__init__()
        self.gamma = nn.Parameter(torch.zeros(1, dim, 1))
        self.beta = nn.Parameter(torch.zeros(1, dim, 1))

    def forward(self, x):
        # norm along spatial dimension then divide
        # that by mean along channel dimension
        Gx = torch.norm(x, p=2, dim=2, keepdim=True)
        Nx = Gx / (Gx.mean(dim=1, keepdim=True) + 1e-5)

        return self.gamma * x * Nx + self.beta + x


class DepthWiseConv1D(nn.Module):
    """
    Depthwise convolution is grouped convolution where the number of groups
    is equal to the number of input channels. This is equivalent to applying
    a different kernel to each channel.
    """

    def __init__(self, in_channels, out_channels, kernel_size):
        super(DepthWiseConv1D, self).__init__()
        self.depthwise = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            groups=in_channels,
            padding="same",
        )

    def forward(self, x):
        return self.depthwise(x)


class PointWiseConv1D(nn.Module):
    """
    Pointwise convolution is a 1x1 convolution that is used to change the number
    of channels in a feature map.
    """

    def __init__(self, in_channels, out_channels):
        super(PointWiseConv1D, self).__init__()
        self.pointwise = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=1,
        )

    def forward(self, x):
        return self.pointwise(x)


class DownsampleConv1D(nn.Module):
    """
    Downsampling applies a convolution with a large stride to reduce the
    spatial dimension of the input.
    """

    def __init__(self, in_channels, out_channels, kernel_size=2):
        super(DownsampleConv1D, self).__init__()

        self.block = nn.Sequential(
            nn.GroupNorm(1, in_channels),
            nn.Conv1d(
                in_channels=in_channels,
                out_channels=out_channels,
                kernel_size=kernel_size,
                stride=kernel_size,
            ),
        )

    def forward(self, x):
        return self.block(x)


class ConvNext1DStem(nn.Module):
    """
    Input stage of cnn that applies a downsample convolution followed by normalization.
    """

    def __init__(self, in_channels=1, out_channels=32, kernel_size=4):
        super(ConvNext1DStem, self).__init__()
        self.stem = nn.Sequential(
            nn.Conv1d(
                in_channels=in_channels,
                out_channels=out_channels,
                kernel_size=kernel_size,
                stride=kernel_size,
            ),
            nn.GroupNorm(1, out_channels),
        )

    def forward(self, x):
        return self.stem(x)


class ConvNext1DBlock(nn.Module):
    """
    Main block of a ConvNext1D network. This block uses the inverted bottleneck
    configuration with skip connection and stochastic depth.
    """

    def __init__(
        self,
        channels,
        kernel_size=7,
        bottleneck_factor=2,
        stride=1,
        stochastic_depth=0.1,
    ):
        super(ConvNext1DBlock, self).__init__()
        self.stochastic_depth = stochastic_depth

        self.block = nn.Sequential(
            DepthWiseConv1D(
                in_channels=channels,
                out_channels=channels,
                kernel_size=kernel_size,
            ),
            nn.GroupNorm(1, channels),
            PointWiseConv1D(
                in_channels=channels,
                out_channels=channels * bottleneck_factor,
            ),
            nn.GELU(),
            GlobalResponseNormalization(channels * bottleneck_factor),
            PointWiseConv1D(
                in_channels=channels * bottleneck_factor,
                out_channels=channels,
            ),
        )

    def forward(self, x):
        """
        Forward pass of ConvNext1D with skip connection and stochastic depth
        (when training). Use the block path with probability {stochastic_depth}
        otherwise use just the identity path.
        """
        # if self.training:
        #     if np.random.rand() <= self.stochastic_depth:
        #         return x
        #     else:
        #         return self.block(x) + x
        # else:
        #     return self.block(x) + x
        return self.block(x) + x


class ConvNext1DStage(nn.Module):
    """
    Stage of ConvNext1D that applies a downsample (if specified)
    followed by n_blocks X ConvNext1DBlock.
    """

    def __init__(
        self,
        in_channels,
        out_channels,
        block_kernel_size=7,
        downsample_kernel_size: int | None = 2,
        n_blocks=3,
        bottleneck_factor=2,
        stochastic_depth=0.1,
    ):
        super(ConvNext1DStage, self).__init__()
        self.blocks = nn.Sequential()
        if downsample_kernel_size:
            self.blocks.append(
                DownsampleConv1D(
                    in_channels=in_channels,
                    out_channels=out_channels,
                    kernel_size=downsample_kernel_size,
                )
            )
        self.blocks.append(
            nn.Sequential(
                *[
                    ConvNext1DBlock(
                        channels=out_channels,
                        kernel_size=block_kernel_size,
                        bottleneck_factor=bottleneck_factor,
                        stochastic_depth=stochastic_depth,
                    )
                    for _ in range(n_blocks - 1)
                ],
                ConvNext1DBlock(
                    channels=out_channels,
                    kernel_size=block_kernel_size,
                    bottleneck_factor=bottleneck_factor,
                    stochastic_depth=0.0,
                ),
            )
        )

    def forward(self, x):
        return self.blocks(x)


class ConvNext1DEncoder(pl.LightningModule):
    def __init__(
        self,
        stem_kernel: int = 4,
        downsample_kernel: int = 2,
        block_kernel: int = 7,
        block_dims: list[int] = [128, 256, 512],
        n_blocks: list[int] = [2, 4, 2],
        stochastic_depths: list[float] = [0.25, 0.5, 0.25],
    ):
        """
        Initialize the ConvNext1D encoder
        Params:
            stem_kernel: kernel size for the stem convolution
            downsample_kernel: kernel size for the downsample convolution, if None, then
                               Then the stage will not downsample (eg in the first stage)
            block_kernel: kernel size for the ConvNext1DBlock depthwise convolution
            block_dims: List of dimensions for each stage.
                        Last element it the encoder dimension.
            n_blocks: List of number of blocks for each stage
            stochastic_depths: List of stochastic depth values for each stage
        """
        super(ConvNext1DEncoder, self).__init__()
        assert (
            len(block_dims) == len(n_blocks) == len(stochastic_depths)
        ), "block_dims, n_blocks, and stochastic_depths must be the same length"

        self.enc_dimension = block_dims[-1]

        # stem and first stage ----------------------------
        self.blocks = nn.Sequential(
            ConvNext1DStem(
                in_channels=1,
                out_channels=block_dims[0],
                kernel_size=stem_kernel,
            ),
        )
        self.blocks.append(
            ConvNext1DStage(
                in_channels=block_dims[0],
                out_channels=block_dims[0],
                block_kernel_size=block_kernel,
                downsample_kernel_size=None,
                n_blocks=n_blocks[0],
                stochastic_depth=stochastic_depths[0],
            )
        )

        # remaining stages --------------------------------
        for i in range(1, len(block_dims)):
            self.blocks.append(
                ConvNext1DStage(
                    in_channels=block_dims[i - 1],
                    out_channels=block_dims[i],
                    block_kernel_size=block_kernel,
                    downsample_kernel_size=downsample_kernel,
                    n_blocks=n_blocks[i],
                    stochastic_depth=stochastic_depths[i],
                )
            )

        # to embedding ------------------------------------
        self.blocks.append(
            nn.Sequential(
                nn.AdaptiveAvgPool1d(1),
                nn.GroupNorm(1, block_dims[-1]),
                nn.Flatten(),
                # TODO test without the linear output layer
                nn.Linear(block_dims[-1], block_dims[-1]),
            )
        )

        self.apply(self._init_weights)

    def forward(self, x):
        return self.blocks(x)
    def predict_step(self, batch, _):
        return self.forward(batch["P"])

    def _init_weights(self, m, w_init=nn.init.trunc_normal_, b_init=nn.init.zeros_):
        """
        The default init for pytorch seems to cause the beginning of training
        to be a bit suboptimal (uses Xavier?).
        """
        if isinstance(m, nn.Conv1d):
            w_init(m.weight)
            if m.bias is not None:
                b_init(m.bias)


# ------------------------------------------------------------------------------
# regular convolutional network blocks
# ------------------------------------------------------------------------------
class Conv1DBlock(nn.Module):
    """Regular convolution with no skip connection"""

    def __init__(
        self,
        in_channels,
        out_channels,
        dropout=0.1,
        kernel_size=3,
        stride=1,
        padding="valid",
    ):
        super(Conv1DBlock, self).__init__()
        self.conv = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            stride=stride,
            padding=padding,
        )
        self.dropout = nn.Dropout(dropout)
        # self.norm = nn.GroupNorm(1, out_channels)
        self.norm = nn.BatchNorm1d(out_channels)
        self.act = nn.GELU()
        self.pool = nn.MaxPool1d(kernel_size=3, stride=1)

    def forward(self, x):
        x = self.conv(x)
        x = self.dropout(x)
        x = self.norm(x)
        x = self.act(x)
        x = self.pool(x)
        return x


class Conv1DEncoder(pl.LightningModule):
    """
    Conv 1d model that takes in a sequence of genotypes (sparse with positional
    channel) and outputs an encoded N-dimensional vector representation.
    """

    def __init__(
        self,
        in_channels=1,
        kernel_size=3,
        stride=1,
        n_layers=4,
        dropout=0.1,
        padding="valid",
        enc_dimension=512,
    ):
        super().__init__()

        self.in_channels = in_channels
        self.enc_dimension = n_layers * 32
        self.kernel_size = kernel_size
        self.stride = stride
        self.n_layers = n_layers
        self.dropout = dropout
        self.padding = padding
        self.enc_dimension = enc_dimension

        self.conv_in = nn.Conv1d(
            in_channels=in_channels,
            out_channels=32,
            kernel_size=kernel_size,
            stride=stride,
            padding=padding,
        )

        self.conv_blocks = nn.ModuleList(
            [
                Conv1DBlock(
                    in_channels=32 * i,
                    out_channels=(i + 1) * 32,
                    dropout=dropout,
                    kernel_size=kernel_size,
                    stride=stride,
                    # padding=padding,
                )
                for i in range(1, n_layers)
            ]
        )

        # used to flatten the output of the conv blocks
        # TODO just use this as the output encoding and see what happens
        self.avg_pool = nn.AdaptiveAvgPool1d(1)

        self.fc = nn.Linear(n_layers * 32, enc_dimension)
        self.fc_dropout = nn.Dropout(dropout)

        self.apply(self._init_weights)

    def forward(self, x):
        x = self.conv_in(x)
        for block in self.conv_blocks:
            x = block(x)
        x = self.avg_pool(x).squeeze(2)
        x = self.fc(x)
        x = self.fc_dropout(x)
        return x

    def predict_step(self, batch, _):
        return self.forward(batch["P"])

    def _init_weights(self, m, w_init=nn.init.trunc_normal_, b_init=nn.init.zeros_):
        """
        The default init for pytorch seems to cause the beginning of training
        to be a bit suboptimal (uses Xavier?).
        """
        if isinstance(m, nn.Conv1d):
            w_init(m.weight)
            if m.bias is not None:
                b_init(m.bias)


encoder_factory = {
    "conv1d": Conv1DEncoder,
    "ConvNext": ConvNext1DEncoder,
}

# ------------------------------------------------------------------------------
# Main siamese network that uses a generic encoder backbone
# ------------------------------------------------------------------------------
class SiameseModule(pl.LightningModule):
    def __init__(
        self,
        *,
        encoder_type,
        encoder_params,
        optimizer,
        optimizer_params,
        scheduler,
        scheduler_params,
        loss_fn,
    ):
        super().__init__()

        self.encoder_type = encoder_type

        self.encoder = encoder_factory[encoder_type](**encoder_params)
        self.enc_dimension = self.encoder.enc_dimension
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.optimizer_params = optimizer_params
        self.scheduler_params = scheduler_params
        self.loss_fn = loss_fn()
        self.save_hyperparameters()

    def forward(self, batch):
        x1 = batch["P1"]
        x2 = batch["P2"]
        d = batch["D"]
        u = self.encoder(x1)
        with torch.no_grad():
            v = self.encoder(x2)

        # see if eps affects numeric instability in mixed precision settings
        # dpred = F.cosine_similarity(u, v, dim=1, eps=1e-9)
        dpred = F.cosine_similarity(u, v, dim=1, eps=1e-6)
        return self.loss_fn(dpred, d)

    def training_step(self, batch, _):
        loss = self.forward(batch)
        self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=False)
        return loss

    def validation_step(self, batch, _):
        loss = self.forward(batch)
        self.log("val_loss", loss, on_step=False, on_epoch=True, prog_bar=False)
        return loss

    def configure_optimizers(self):
        optimizer = self.optimizer(self.parameters(), **self.optimizer_params)
        scheduler = self.scheduler(optimizer, **self.scheduler_params)
        return {
            "optimizer": optimizer,
            "lr_scheduler": scheduler,
            "monitor": "val_loss",
        }
