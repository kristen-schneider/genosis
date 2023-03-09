import argparse
import os
import random
from datetime import datetime
from functools import partial
from glob import glob

import numpy as np
import pytorch_lightning as pl
import torch
from data_utils.gt_datasets import (GTDataset, GTRandomAugmentDataset,
                                    pad_data, random_augment_pairs,
                                    train_val_split)
from models.encoder import Conv1DEncoder, SiameseModule, SimSiamModule
from pytorch_lightning.callbacks import (EarlyStopping, ModelCheckpoint,
                                         StochasticWeightAveraging)
from pytorch_lightning.loggers.wandb import WandbLogger
from torch import nn, optim
from torch.utils import data
from torchmetrics import MeanSquaredLogError
from torchvision import ops


def compute_steps_per_epoch(train_dataset, batch_size):
    """
    Compute the total number of optimizer steps for the learning rate scheduler.
    """
    n_batches = int(np.ceil(len(train_dataset) / batch_size))
    return n_batches


def get_dataloaders(args):
    # make main dataset
    train_ds = GTDataset(
        pos_files=(args.P1_train, args.P2_train),
        dist_file=args.D_train,
    )
    val_ds = GTDataset(
        pos_files=(args.P1_val, args.P2_val),
        dist_file=args.D_val,
    )

    collate_fn = partial(pad_data, model_type=args.model_type)

    train_dataloader = data.DataLoader(
        train_ds,
        batch_size=args.batch_size,
        collate_fn=collate_fn,
        shuffle=True,
        num_workers=args.n_workers,
        pin_memory=True,
        drop_last=False,
        prefetch_factor=2,
    )
    val_dataloader = data.DataLoader(
        val_ds,
        batch_size=args.batch_size,
        collate_fn=collate_fn,
        shuffle=False,
        num_workers=args.n_workers,
        pin_memory=True,
        drop_last=False,
        prefetch_factor=2,
    )

    return {
        "train_dataloader": train_dataloader,
        "val_dataloader": val_dataloader,
    }


def siamese(args):

    if args.loss_fn == "mse":
        loss_fn = nn.MSELoss
    elif args.loss_fn == "msle":
        loss_fn = MeanSquaredLogError
    elif args.loss_fn == "huber":
        loss_fn = nn.SmoothL1Loss
    else:
        raise ValueError(f"Loss function {args.loss_fn} not supported.")

    data = get_dataloaders(args)

    if args.train_method == "sim_siam":
        cls = SimSiamModule
    else:
        cls = SiameseModule

    siamese_model = cls(
        encoder_type=args.model_type,
        encoder_params={  # TODO add more for transformer if needed
            "n_layers": args.n_layers,
            "dropout": args.dropout,
            "kernel_size": args.kernel_size,
            "stride": args.stride,
            "padding": args.padding,
            # "enc_dimension": args.enc_dimension,
        },
        lr=args.lr,
        optimizer=optim.AdamW,
        optimizer_params={
            # "lr": args.lr,
            "weight_decay": args.weight_decay,
            "eps": 1e-4, # TODO see if we need to lower this even further
        },
        scheduler=optim.lr_scheduler.CosineAnnealingWarmRestarts,
        scheduler_params={
            "T_0": 2
            * compute_steps_per_epoch(
                train_dataset=data["train_dataloader"].dataset,
                batch_size=args.batch_size,
            ),
            "T_mult": 1,
            "eta_min": 1e-6,
            "verbose": True,
        },
        loss_fn=loss_fn,
    )

    callbacks = [
        StochasticWeightAveraging(swa_lrs=args.lr),
        EarlyStopping(
            monitor="val_loss",
            patience=args.early_stop_patience,
            verbose=True,
            mode="min",
        ),
        # save the last 10 checkpoints
        ModelCheckpoint(
            monitor="val_loss",
            dirpath=f"{args.outdir}/{args.prefix}.checkpoints",
            filename=f"{args.train_method}-{args.model_type}-{{epoch:03d}}-{{val_loss:.5f}}",
            save_top_k=10,
            save_last=True,
            mode="min",
        ),
    ]

    logger = WandbLogger(
        project="gt-similarity-search",
        name=f"conv1d-siamese-{args.prefix}-{datetime.now().strftime('%Y-%m-%d-%H-%M')}",
        log_model="all",
    )

    trainer = pl.Trainer(
        precision=16, # TODO add to config
        accumulate_grad_batches=args.grad_accum,
        gradient_clip_val=0.5,
        gradient_clip_algorithm="norm",
        max_epochs=args.n_epochs,
        val_check_interval=0.5,
        callbacks=callbacks,
        logger=logger,
        accelerator="gpu",
        devices=1,
        fast_dev_run=args.fast_dev_run,  # run 1 batch train/val to see if things are working
        enable_progress_bar=False,
    )

    logger.watch(siamese_model)

    torch.set_float32_matmul_precision("high") # I think medium|high|very_high are the options
    trainer.fit(
        model=siamese_model,
        train_dataloaders=data["train_dataloader"],
        val_dataloaders=data["val_dataloader"],
    )


def transformer_paired_lm_pretrain(args):
    """
    Transformer with paired example input and language model objective
    """
    raise NotImplementedError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fast_dev_run",
        action="store_true",
        help="run 1 batch train/val to see if things are working",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="directory to save checkpoints to",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="",
        help="prefix to add to run name and checkpoint directory name",
    )
    parser.add_argument(
        "--P1_train",
        type=str,
        required=True,
        help="Path to training P1 memmap",
    )
    parser.add_argument(
        "--P2_train",
        type=str,
        required=True,
        help="Path to training P2 memmap",
    )
    parser.add_argument(
        "--P1_val",
        type=str,
        required=True,
        help="Path to validation P1 memmap",
    )
    parser.add_argument(
        "--P2_val",
        type=str,
        required=True,
        help="Path to validation P2 memmap",
    )
    parser.add_argument(
        "--D_train",
        type=str,
        required=True,
        help="Path to training distance file",
    )
    parser.add_argument(
        "--D_val",
        type=str,
        required=True,
        help="Path to validation distance file",
    )
    parser.add_argument(
        "--train_method",
        type=str,
        default="sim_siam",
        help="training method",
        choices=[
            "sim_siam",
            "siamese",
            "transformer_paired_lm_pretrain",
        ],
    )
    parser.add_argument(
        "--model_type",
        type=str,
        default="conv1d_siamese",
        help="model type to train",
        choices=[
            "conv1d",
            "transformer",
        ],
    )
    parser.add_argument(
        "--loss_fn",
        type=str,
        default="mse",
        help="loss function to use",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=32,
        help="Batch size for training",
    )
    parser.add_argument(
        "--grad_accum",
        type=int,
        default=1,
        help="Number of gradient accumulation steps",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=4,
        help="Number of workers for dataloader",
    )
    parser.add_argument(
        "--n_epochs",
        type=int,
        default=10,
        help="Number of epochs to train for",
    )
    parser.add_argument(
        "--early_stop_patience",
        type=int,
        default=5,
        help="Number of validation checks to wait before early stopping",
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=1e-3,
        help="Base learning rate for optimizer",
    )
    parser.add_argument(
        "--weight_decay",
        type=float,
        default=1e-3,
        help="Weight decay for optimizer (if applicable)",
    )
    parser.add_argument(
        "--n_layers",
        type=int,
        default=4,
        help="Number of layers in the encoder",
    )
    parser.add_argument(
        "--dropout",
        type=float,
        default=0.1,
        help="Dropout rate for encoder (if applicable)",
    )
    parser.add_argument(
        "--kernel_size",
        type=int,
        default=3,
        help="Kernel size (for CNNs only)",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Stride (for CNNs only)",
    )
    parser.add_argument(
        "--padding",
        type=str,
        default="valid",
        help="Padding (for CNNs only)",
    )
    args = parser.parse_args()

    if args.model_type == "conv1d" or args.model_type == "transformer_siamese":
        siamese(args)
    elif args.model_type == "transformer_paired_lm_pretrain":
        transformer_paired_lm_pretrain(args)
    else:
        raise ValueError("Invalid model type")


if __name__ == "__main__":
    main()
