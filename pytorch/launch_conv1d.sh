
FAST_DEV_RUN=$([[ $2 == "fast_dev_run" ]] && echo "--fast_dev_run" || echo "")

python train_model.py \
    --prefix $1 \
    --P1 "/data/P1" \
    --P2 "/data/P2" \
    --D "/data/precomputed_distances/D.txt" \
    --train_method siamese \
    --model_type conv1d \
    --batch_size 64 \
    --grad_accum 1 \
    --n_workers 4 \
    --n_epochs 100 \
    --early_stop_patience 50 \
    --lr 0.001 \
    --weight_decay 0.001 \
    --n_layers 3 \
    --dropout 0.0 \
    --kernel_size 3 \
    --stride 1 \
    --padding same \
    $FAST_DEV_RUN
