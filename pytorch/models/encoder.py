import numpy as np
import pytorch_lightning as pl
import torch
import torch.nn.functional as F
from torch import nn
from transformers.models.longformer.configuration_longformer import \
    LongformerConfig
from transformers.models.longformer.modeling_longformer import LongformerModel


class ResidualBlock(nn.Module):
    """conv1d block with residual connection"""

    def __init__(
        self,
        in_channels,
        out_channels,
        dropout=0.1,
        kernel_size=3,
        stride=1,
    ):
        super(ResidualBlock, self).__init__()

        self.in_channels = in_channels
        self.out_channels = out_channels

        self.conv1 = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            stride=stride,
            padding=1,
        )
        self.conv2 = nn.Conv1d(
            in_channels=out_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            stride=stride,
            padding=1,
        )
        # allows residual connection to be applied when in_channels != out_channels
        self.downsample = nn.Sequential(
            nn.Conv1d(
                in_channels=in_channels,
                out_channels=out_channels,
                kernel_size=1,
                stride=stride,
            ),
            nn.BatchNorm1d(out_channels),
            # nn.GroupNorm(1, out_channels),
        )
        self.dropout = nn.Dropout(dropout)
        # self.norm = nn.GroupNorm(1, out_channels)
        self.norm = nn.BatchNorm1d(out_channels)
        self.act1 = nn.GELU()
        self.act2 = nn.GELU()
        # self.pool = nn.MaxPool1d(kernel_size=3, stride=1)

    def forward(self, x):
        identity = x
        out = self.conv1(x)
        out = self.dropout(out)
        out = self.norm(out)
        out = self.act1(out)
        out = self.conv2(out)
        out = self.dropout(out)
        out = self.norm(out)
        if self.in_channels != self.out_channels:
            identity = self.downsample(identity)
        out += identity
        out = self.act2(out)
        # out = self.pool(out)
        return out


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

    def forward(self, x):
        x = self.conv_in(x)
        for block in self.conv_blocks:
            x = block(x)
        x = self.avg_pool(x).squeeze(2)
        x = self.fc(x)
        x = self.fc_dropout(x)
        return x
        # return F.normalize(x, dim=1)

    def predict_step(self, batch, _):
        return self.forward(batch["P"])


class SimSiamModule(pl.LightningModule):
    def __init__(
        self,
        *,
        encoder_type,
        encoder_params,
        lr,
        optimizer,
        optimizer_params,
        scheduler,
        scheduler_params,
        loss_fn=None,  # not used
    ):
        super().__init__()

        self.encoder_type = encoder_type
        self.encoder = Conv1DEncoder(**encoder_params)
        self.enc_dimension = self.encoder.enc_dimension
        self.projection_head = nn.Sequential(
            nn.Linear(self.enc_dimension , self.enc_dimension // 2),
            nn.BatchNorm1d(self.enc_dimension // 2),
            nn.ReLU(),
            nn.Linear(self.enc_dimension // 2, self.enc_dimension),
        )
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.optimizer_params = optimizer_params
        self.scheduler_params = scheduler_params

        # 2-layer projection head
        # TODO todo make the input and output satisfy the encoder params in general

        if lr:
            # precedence to lr passed in
            self.optimizer_params["lr"] = lr
            self.learning_rate = lr
        self.save_hyperparameters()

    def forward(self, x):
        p1 = self.encoder(x["P1"])
        z1 = self.projection_head(p1)

        p2 = self.encoder(x["P2"])
        z2 = self.projection_head(p2)
        return cos_sim_loss(p1.detach(), z2)  # + cos_sim_loss(p2.detach(), z1))

    def training_step(self, batch, _):
        loss = self.forward(batch)
        self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=True)
        return loss

    def validation_step(self, batch, _):
        loss = self.forward(batch)
        self.log("val_loss", loss, prog_bar=True)
        return loss

    def configure_optimizers(self):
        optimizer = self.optimizer(self.parameters(), **self.optimizer_params)
        scheduler = self.scheduler(optimizer, **self.scheduler_params)
        # return [optimizer], [scheduler]
        return {
            "optimizer": optimizer,
            "lr_scheduler": scheduler,
            "monitor": "val_loss",
        }


# TODO change the args so that we can load from checkpoint without specifying it
class SiameseModule(pl.LightningModule):
    def __init__(
        self,
        *,
        encoder_type,
        encoder_params,
        lr,
        optimizer,
        optimizer_params,
        scheduler,
        scheduler_params,
        loss_fn,
    ):
        super().__init__()

        self.encoder_type = encoder_type

        if self.encoder_type == "conv1d":
            self.encoder = Conv1DEncoder(**encoder_params)
            self.enc_dimension = self.encoder.enc_dimension
        else:
            raise ValueError(f"Model type {self.encoder_type} not supported.")



        self.optimizer = optimizer
        self.scheduler = scheduler
        self.optimizer_params = optimizer_params
        self.scheduler_params = scheduler_params
        self.loss_fn = loss_fn()
        if lr:
            # precedence to lr passed in
            self.optimizer_params["lr"] = lr
            self.learning_rate = lr
        self.save_hyperparameters()

    def forward(self, batch):
        x1 = batch["P1"]
        x2 = batch["P2"]
        d = batch["D"]
        # TODO clean this up
        u = self.encoder(x1)
        with torch.no_grad():
            v = self.encoder(x2)

        dpred = F.cosine_similarity(u, v, dim=1, eps=1e-9)
        return self.loss_fn(dpred, d)

    def training_step(self, batch, _):
        loss = self.forward(batch)
        self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=True)
        return loss

    def validation_step(self, batch, _):
        x1 = batch["P1"]
        x2 = batch["P2"]
        d = batch["D"]
        with torch.no_grad():
            u = self.encoder(x1)
            v = self.encoder(x2)
            dpred = F.cosine_similarity(u, v, dim=1, eps=1e-9)
            loss = F.mse_loss(dpred, d)
        self.log("val_loss", loss, prog_bar=True)
        return loss

    def configure_optimizers(self):
        optimizer = self.optimizer(self.parameters(), **self.optimizer_params)
        scheduler = self.scheduler(optimizer, **self.scheduler_params)
        # return [optimizer], [scheduler]
        return {
            "optimizer": optimizer,
            "lr_scheduler": scheduler,
            "monitor": "val_loss",
        }
