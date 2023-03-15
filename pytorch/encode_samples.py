import argparse
import os
from functools import partial
from pprint import pprint

import pytorch_lightning as pl
import torch
from data_utils.gt_datasets import (GTInferenceDataset, GTInferenceWriter,
                                    pad_data)
from models.encoder import SiameseModule
from torch.utils.data import DataLoader


def load_siamese_encoder(path: str) -> pl.LightningModule:
    # model = SiameseModule.load_from_checkpoint(path)
    model = SiameseModule.load_from_checkpoint(path)
    encoder = model.encoder
    return encoder


def encode_samples(
    *,
    encoder: pl.LightningModule,
    batch_size: int,
    output: str,
    files: list[str],
    gpu: bool = False,
    num_workers: int = 0,
):
    """
    Encode a dataset of samples using the provided encoder.
    :param encoder: The encoder model for inference
    :param batch_size: The batch size for inference
    :param output: The output file to write the encoded samples to
    :param files: The files to encode
    :param gpu: Whether to use the GPU
    :param num_workers: The number of workers to use for data loading
    """
    pprint(files)

    dataset = GTInferenceDataset(files=files)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=True,
        collate_fn=partial(pad_data, model_type="conv1d_inference"),
    )
    writer = GTInferenceWriter(
        output=output,
        write_interval="batch",
    )
    trainer = pl.Trainer(
        inference_mode=True,
        callbacks=[writer],
        accelerator="cuda" if gpu else "cpu",
        devices=1 if gpu else None,
    )
    trainer.predict(encoder, dataloader)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--encoder",
        type=str,
        required=True,
        help="Path to the encoder checkpoint",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="output file",
    )
    parser.add_argument(
        "--files",
        type=str,
        nargs="+",
        required=True,
        help="Files to encode",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size for inference",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        default=False,
        help="Use GPU",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help="Number of workers for data loading.",
    )
    args = parser.parse_args()
    encoder = load_siamese_encoder(args.encoder)
    encode_samples(
        encoder=encoder,
        batch_size=args.batch_size,
        output=args.output,
        files=args.files,
        gpu=args.gpu,
        num_workers=args.num_workers,
    )
