python create_vectors.py \
    --sample-list ~/precision-medicine/data/samples/chr8-30x/testing.samples \
    --model-path base_model.h5 \
    --output-path ~/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.0.embeddings \
    --batch-size 32 \
    --sample_id_filename ~/precision-medicine/data/samples/chr8-30x/all.samples \
    --genotype_filename ~/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.0.encoded
