"""
Contains various datasets/dataloaders for use with GT similarity models.
"""
from functools import partial
from glob import glob
from itertools import chain
from pprint import pprint
from typing import Sequence

import numpy as np
import torch
import torch.utils.data as data
from mmap_ninja.ragged import RaggedMmap
from pytorch_lightning.callbacks import BasePredictionWriter
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, Dataset, IterableDataset, Subset

# (nsamples in 1000 genomes) * (2 for each allele)
NSAMPLES = 6404

# TODO add a way to do this directly from a vcf?
class GTInferenceDataset(IterableDataset):
    def __init__(
        self,
        files: list[str],  # list of positional encodings
    ):
        super().__init__()
        self.files = files
        self.n_files = len(files)

    def _make_iter(self, idx):
        """
        Make an iterator over one file.
        """
        with open(self.files[idx], "r") as f:
            # filename format 1KG.data.seg.NNN.pos_encoded
            # TODO later need to make this more general
            #segment = self.files[idx].split(".")[-2]
            segment = self.files[idx].split('.')[-2].replace('segment','')
            for line in f:
                # sample name, then 1D array of positional encodings
                A = line.rstrip().split()
                yield {
                    "segment": segment,
                    "sample": A[0],
                    "P": np.array(A[1:], dtype=np.float32) - float(segment),
                }

    def __iter__(self):
        workers_info = data.get_worker_info()
        if workers_info is None:
            # single process
            iter_start = 0
            iter_end = self.n_files
        else:
            # multiprocessing: compute division of labor
            per_worker = int(np.ceil(self.n_files / workers_info.num_workers))
            worker_id = workers_info.id
            iter_start = worker_id * per_worker
            iter_end = min(iter_start + per_worker, self.n_files)

        return chain.from_iterable(
            self._make_iter(i) for i in range(iter_start, iter_end)
        )


class GTInferenceWriter(BasePredictionWriter):
    """
    Write the results of inference to a formatted file.
    Format:
    $SAMPLE $SEGMENT $EMBEDDING (tab separated, newline terminated)
    """

    def __init__(self, *, output, write_interval):
        super().__init__(write_interval)
        self.output = output
        # if the file already exists, delete it
        with open(f"{self.output}", "w") as f:
            f.write("")

    def write_on_batch_end(
        self,
        trainer,
        pl_module,
        outputs,
        batch_indices,
        batch,
        batch_idx,
        dataloader_idx,
    ):
        if trainer.num_devices > 0:
            outputs = outputs.cpu()

        with open(f"{self.output}", "a") as f:
            for sample, segment, embedding in zip(
                batch["sample"], batch["segment"], outputs
            ):
                f.write(
                    f"{sample}\t{segment}\t{' '.join(map(str, embedding.numpy()))}\n"
                )


class GTRandomAugmentDataset(Dataset):
    def __init__(self, pos_file: str):
        self.P = RaggedMmap(pos_file)
        self.n_data = len(self.P)

    def __getitem__(self, idx):
        return self.P[idx]

    def __len__(self):
        return self.n_data


def random_augment_pairs(batch, p=0.1):
    """
    - data: list of tensors of positional encodings.
    - p: probability of omitting a position

    Takes a batch of positional encodings and returns a batch of
    pairs: (P1, P2), where P1 and P2 are randomly augmented versions
    of the same positional encoding.

    The random augmentation will be created by randomly ommitting
    a certain number of positions with probability p.

    Additionally, all P1/P2 tensors will be padded to the same length
    of their respective batches.
    """
    # copy over batch to arrays
    P = [np.array(x) for x in batch]
    P_noise1 = [x + np.random.normal(0, 0.001, size=x.shape) for x in P]
    P_noise2 = [x + np.random.normal(0, 0.001, size=x.shape) for x in P]

    # randomly omit positions of P with probability p
    P1_aug = [torch.Tensor(x[np.random.rand(len(x)) > p]) for x in P_noise1]
    P2_aug = [torch.Tensor(x[np.random.rand(len(x)) > p]) for x in P_noise2]


    # Pad all P1/P2 tensors to the same length (separately).
    # Randomly choose which tensor is the augmented one.
    P1 = pad_sequence(P1_aug, batch_first=True).unsqueeze(1)
    P2 = pad_sequence(P2_aug, batch_first=True).unsqueeze(1)
    return {"P1": P1, "P2": P2}


class GTDataset(Dataset):
    """
    Load the gentoypes and cM positions and create batches of tensors
    """

    def __init__(
        self,
        *,
        pos_files: Sequence[str],
        dist_file: str,
    ):
        self.dist_file = dist_file

        self.P1 = RaggedMmap(pos_files[0])
        self.P2 = RaggedMmap(pos_files[1])
        self.D = np.array([float(x.rstrip()) for x in open(self.dist_file, "r")])
        self.n_data = len(self.D)
        assert len(self.P1) == len(self.P2) == self.n_data

    def __getitem__(self, idx):
        return {"P1": self.P1[idx], "P2": self.P2[idx], "D": self.D[idx]}

    def __len__(self):
        return self.n_data


# TODO maybe not the best way to do this
# could just separate out the functions
def pad_data(data, model_type="conv1d"):
    # TODO this is fragile
    if model_type == "conv1d" or model_type == "ConvNext":
        # variable length tensors
        P1 = [torch.tensor(x["P1"], dtype=torch.float32) for x in data]
        P2 = [torch.tensor(x["P2"], dtype=torch.float32) for x in data]

        # scalars
        dists = torch.tensor([x["D"] for x in data], dtype=torch.float32)

        P1 = pad_sequence(P1, batch_first=True).unsqueeze(1)
        P2 = pad_sequence(P2, batch_first=True).unsqueeze(1)
        return {"P1": P1, "P2": P2, "D": dists}

    elif model_type == "conv1d_inference":
        P = pad_sequence(
            [torch.tensor(x["P"], dtype=torch.float32) for x in data], batch_first=True
        ).unsqueeze(1)
        samples = [x["sample"] for x in data]
        segments = [x["segment"] for x in data]

        return {"P": P, "sample": samples, "segment": segments}

    elif model_type == "transformer_siamese":
        # variable length tensors
        P1 = [torch.tensor(x["P1"], dtype=torch.float32) for x in data]
        P2 = [torch.tensor(x["P2"], dtype=torch.float32) for x in data]

        # scalars
        dists = torch.tensor([x["D"] for x in data], dtype=torch.float32)

        P1 = pad_sequence(P1, batch_first=True).unsqueeze(1)
        P2 = pad_sequence(P2, batch_first=True).unsqueeze(1)

        attention_mask1 = torch.where(P1 == 0, 0, 1).squeeze(2)
        attention_mask2 = torch.where(P2 == 0, 0, 1).squeeze(2)
        return {
            "P1": P1,
            "P2": P2,
            "attention_mask1": attention_mask1,
            "attention_mask2": attention_mask2,
            "D": dists,
        }

    else:
        raise ValueError(f"Unknown model type: {model_type}")


def train_val_split(
    dataset: Dataset,
    train_segments: Sequence[int],
    val_segments: Sequence[int],
    nsamples: int = NSAMPLES,
):
    """
    Subset the dataset into train and validation sets with indices provided by the user.
    """
    assert len(dataset) % nsamples == 0
    train_dataset = Subset(dataset, train_segments)
    val_dataset = Subset(dataset, val_segments)
    return train_dataset, val_dataset


# ==============================================================================
# Sanity Check
# ==============================================================================
if __name__ == "__main__":
    # files = glob("/data/segments/*.pos_encoded")[:10]

    # files = [
    #     "/data/segments/1KG.data.seg.90.pos_encoded",
    #     "/data/segments/1KG.data.seg.91.pos_encoded",
    #     "/data/segments/1KG.data.seg.92.pos_encoded",
    #     "/data/segments/1KG.data.seg.93.pos_encoded",
    #     "/data/segments/1KG.data.seg.94.pos_encoded",
    #     "/data/segments/1KG.data.seg.95.pos_encoded",
    #     "/data/segments/1KG.data.seg.96.pos_encoded",
    #     "/data/segments/1KG.data.seg.97.pos_encoded",
    #     "/data/segments/1KG.data.seg.98.pos_encoded",
    #     "/data/segments/1KG.data.seg.99.pos_encoded",
    #     "/data/segments/1KG.data.seg.9.pos_encoded",
    # ]
    # infer_dataset = GTInferenceDataset(files=files)

    dataset = GTDataset(
        pos_files=("/data/P1", "/data/P2"),
        dist_file="/data/precomputed_distances/D.txt",
        model_type="transformer_siamese",
    )

    print(max([len(x["P1"]) for x in dataset]))
