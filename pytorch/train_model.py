import argparse
import os
import random
from datetime import datetime
from glob import glob

import numpy as np
import pytorch_lightning as pl
from data_utils.gt_datasets import GTDataset, pad_data, train_val_split
from models.encoder import Conv1DEncoder, SiameseModule
from pytorch_lightning.callbacks import (EarlyStopping, ModelCheckpoint,
                                         StochasticWeightAveraging)
from pytorch_lightning.loggers.wandb import WandbLogger
from torch import nn, optim
from torch.utils import data


def compute_steps_per_epoch(train_dataset, batch_size):
    """
    Compute the total number of optimizer steps for the learning rate scheduler.
    """
    n_batches = int(np.ceil(len(train_dataset) / batch_size))
    return n_batches


def get_dataloaders(args):
    # make main dataset
    dataset = GTDataset(
        pos_files=(args.P1, args.P2),
        dist_file=args.D,
        model_type=args.model_type,
    )

    N = len(dataset)
    train_ds, val_ds = train_val_split(
        train_segments=range(N - 10),
        val_segments=range(N - 10, N - 5),
        dataset=dataset,
    )

    train_dataloader = data.DataLoader(
        train_ds,
        batch_size=args.batch_size,
        collate_fn=pad_data,
        shuffle=True,
        num_workers=args.n_workers,
        pin_memory=True,
        drop_last=False,
        prefetch_factor=2,
    )
    val_dataloader = data.DataLoader(
        val_ds,
        batch_size=args.batch_size,
        collate_fn=pad_data,
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


def conv1d_siamese(args):

    data = get_dataloaders(args)

    siamese_model = SiameseModule(
        encoder_type=args.model_type,
        encoder_params={
            "n_layers": args.n_layers,
            "dropout": args.dropout,
            "kernel_size": args.kernel_size,
            "stride": args.stride,
            "padding": args.padding,
        },
        lr=args.lr,
        optimizer=optim.AdamW,
        optimizer_params={
            # "lr": args.lr,
            "weight_decay": args.weight_decay,
        },
        # scheduler=optim.lr_scheduler.ReduceLROnPlateau,
        # scheduler_params={
        #     "mode": "min",
        #     "factor": 0.1,
        #     "patience": 2,
        #     "verbose": True,
        # },
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
        loss_fn=nn.MSELoss,
    )

    callbacks = [
        StochasticWeightAveraging(swa_lrs=args.lr),
        EarlyStopping(
            monitor="val_loss",
            patience=args.early_stop_patience,
            verbose=True,
            mode="min",
        ),
        ModelCheckpoint(
            monitor="val_loss",
            dirpath=f"{args.prefix}.checkpoints",
            filename="siamese-{epoch:02d}-{val_loss:.2f}",
            save_top_k=3,
            mode="min",
        ),
    ]

    logger = WandbLogger(
        project="gt-similarity-search",
        name=f"conv1d-siamese-{args.prefix}-{datetime.now().strftime('%Y-%m-%d-%H-%M')}",
        log_model="all",
    )

    trainer = pl.Trainer(
        accumulate_grad_batches=args.grad_accum,
        gradient_clip_val=0.5,
        gradient_clip_algorithm="norm",
        max_epochs=args.n_epochs,
        val_check_interval=0.25,
        callbacks=callbacks,
        logger=logger,
        accelerator="gpu",
        devices=1,
        fast_dev_run=args.fast_dev_run,  # run 1 batch train/val to see if things are working
    )

    logger.watch(siamese_model)

    trainer.fit(
        model=siamese_model,
        train_dataloaders=data["train_dataloader"],
        val_dataloaders=data["val_dataloader"],
    )


def transformer_siamese(args):
    """
    Use a transformer as the encoder in a siamese network.
    """
    raise NotImplementedError


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
        "--prefix",
        type=str,
        default="",
        help="prefix to add to run name and checkpoint directory name",
    )
    parser.add_argument(
        "--P1",
        type=str,
        required=True,
        help="Path to P1 memmap",
    )
    parser.add_argument(
        "--P2",
        type=str,
        required=True,
        help="Path to P2 memmap",
    )
    parser.add_argument(
        "--D",
        type=str,
        required=True,
        help="Path to distance file",
    )
    parser.add_argument(
        "--val_segments",
        type=int,
        default=5,
        help="number of validation segments to use",
    )
    parser.add_argument(
        "--model_type",
        type=str,
        default="conv1d_siamese",
        help="model type to train",
        choices=[
            "conv1d_siamese",
            "transformer_siamese",
            "transformer_paired_lm_pretrain",
        ],
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

    if args.model_type == "conv1d_siamese":
        conv1d_siamese(args)
    elif args.model_type == "transformer_siamese":
        transformer_siamese(args)
    elif args.model_type == "transformer_paired_lm_pretrain":
        transformer_paired_lm_pretrain(args)
    else:
        raise ValueError("Invalid model type")


if __name__ == "__main__":
    main()
