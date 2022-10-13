PYTHON_DIR="/home/sdp/precision-medicine/python/scripts/faiss/"
FAISS_ENC="/home/sdp/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.86.enc.faiss"
FAISS_EMB="/home/sdp/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.86.emb.faiss"
PLINK="/home/sdp/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.86.genome"
EUC_DIST="/home/sdp/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.86.test.edist"
TEST="/home/sdp/precision-medicine/data/samples/chr8-30x/testing.samples"
TRAIN="/home/sdp/precision-medicine/data/samples/chr8-30x/training.samples"

python $PYTHON_DIR"performance.py" \
	--faiss_enc $FAISS_ENC\
	--faiss_emb $FAISS_EMB \
	--plink $PLINK \
	--ed $EUC_DIST \
	--test $TEST \
	--train $TRAIN
