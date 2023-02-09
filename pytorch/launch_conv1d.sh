# usage: train_model.py [-h] [--fast_dev_run] [--prefix PREFIX] --P1 P1 --P2 P2 --D D
#                       [--val_segments VAL_SEGMENTS]
#                       [--model_type {conv1d_siamese,transformer_siamese,transformer_paired_lm_pretrain}]
#                       [--batch_size BATCH_SIZE] [--grad_accum GRAD_ACCUM]
#                       [--n_workers N_WORKERS] [--n_epochs N_EPOCHS]
#                       [--early_stop_patience EARLY_STOP_PATIENCE] [--lr LR]
#                       [--weight_decay WEIGHT_DECAY] [--n_layers N_LAYERS] [--dropout DROPOUT]
#                       [--kernel_size KERNEL_SIZE] [--stride STRIDE] [--padding PADDING]

# options:
#   -h, --help            show this help message and exit
#   --fast_dev_run        run 1 batch train/val to see if things are working
#   --prefix PREFIX       prefix to add to run name and checkpoint directory name
#   --P1 P1               Path to P1 memmap
#   --P2 P2               Path to P2 memmap
#   --D D                 Path to distance file
#   --val_segments VAL_SEGMENTS
#                         number of validation segments to use
#   --model_type {conv1d_siamese,transformer_siamese,transformer_paired_lm_pretrain}
#                         model type to train
#   --batch_size BATCH_SIZE
#                         Batch size for training
#   --grad_accum GRAD_ACCUM
#                         Number of gradient accumulation steps
#   --n_workers N_WORKERS
#                         Number of workers for dataloader
#   --n_epochs N_EPOCHS   Number of epochs to train for
#   --early_stop_patience EARLY_STOP_PATIENCE
#                         Number of validation checks to wait before early stopping
#   --lr LR               Base learning rate for optimizer
#   --weight_decay WEIGHT_DECAY
#                         Weight decay for optimizer (if applicable)
#   --n_layers N_LAYERS   Number of layers in the encoder
#   --dropout DROPOUT     Dropout rate for encoder (if applicable)
#   --kernel_size KERNEL_SIZE
#                         Kernel size (for CNNs only)
#   --stride STRIDE       Stride (for CNNs only)
#   --padding PADDING     Padding (for CNNs only)


FAST_DEV_RUN=$([[ $2 == "fast_dev_run" ]] && echo "--fast_dev_run" || echo "")

python train_model.py \
    --prefix $1 \
    --P1 "/data/P1" \
    --P2 "/data/P2" \
    --D "/data/precomputed_distances/D.txt" \
    --val_segments 5 \
    --model_type conv1d_siamese \
    --batch_size 64 \
    --grad_accum 1 \
    --n_workers 4 \
    --n_epochs 100 \
    --early_stop_patience 20 \
    --lr 0.001 \
    --weight_decay 0.01 \
    --n_layers 3 \
    --dropout 0.0 \
    --kernel_size 3 \
    --stride 1 \
    --padding same \
    $FAST_DEV_RUN
