import tensorflow as tf
from datasets import PairsDataset, partition_samples
from models import SiameseModel, build_siamese_network, resnet_model

# load the dataset
# TODO split the data into train and test at the file level
# - do the split by sample name
# - take 90% of samples for training and 5%/5% for val/test
# - use grep -v to hold out samples from the training/val/test sets

if __name__ == "__main__":

    train_samples, val_samples, test_samples = partition_samples(
        train_ratio=0.9, sample_file="/home/murad/data/toy_model_data/ALL.sampleIDs"
    )
    # log the train/val/test samples
    with open("train_samples.txt", "w") as f:
        f.write("\n".join(train_samples))
    with open("val_samples.txt", "w") as f:
        f.write("\n".join(val_samples))
    with open("test_samples.txt", "w") as f:
        f.write("\n".join(test_samples))

    # TODO test to make sure there is no train/test leakage
    # TODO also the filtering is slow.  Not  a big deal, but I would
    # like to be able to do it faster.
    train_data = PairsDataset(
        keep_samples=train_samples,
        sample_id_filename="/home/murad/data/toy_model_data/ALL.sampleIDs",
        sample_pair_filename="/home/murad/data/toy_model_data/chr0.seg.0.cnn",
        genotype_filename="/home/murad/data/toy_model_data/chr0.seg.0.encoding",
        shuffle=True,
        batch_size=128,
    )
    val_data = PairsDataset(
        keep_samples=val_samples,
        sample_id_filename="/home/murad/data/toy_model_data/ALL.sampleIDs",
        sample_pair_filename="/home/murad/data/toy_model_data/chr0.seg.0.cnn",
        genotype_filename="/home/murad/data/toy_model_data/chr0.seg.0.encoding",
        shuffle=True,
        batch_size=128,
    )

    # get base model
    # TODO just wrap this into the siamze model class?
    base_model = resnet_model(shape=(train_data.num_variants, 1))
    # sanity check to compare with trained model
    base_model.save("/home/murad/data/toy_model_data/base_model_untrained.h5")

    siamese_net = build_siamese_network(
        base_model=base_model,
        vector_size=train_data.num_variants,
    )
    training_model = SiameseModel(siamese_net)
    training_model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
        run_eagerly=True,
    )

    callbacks = [
        tf.keras.callbacks.EarlyStopping(
            monitor="val_loss", patience=2, verbose=1
        ),
        # TODO DOESN'T WORK
        # find another callback that saves just the 'base_model'.
        # for now just save at the end
        # tf.keras.callbacks.ModelCheckpoint(
        #     filepath="/home/murad/data/toy_model_data/checkpoints/checkpoint.tf",
        #     monitor="val_loss",
        #     save_best_only=True,
        #     verbose=1,
        #     format="tf"
        # ),
    ]

    training_model.fit(
        train_data.ds,
        epochs=10,
        steps_per_epoch=train_data.num_pairs,
        validation_data=val_data.ds,
        validation_steps=val_data.num_pairs,
        callbacks=callbacks,
    )
    base_model.save("/home/murad/data/toy_model_data/base_model.h5")
