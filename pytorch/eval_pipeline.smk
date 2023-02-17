import os
import sys
from types import SimpleNamespace

"""
# DONE make sure that sample ids are sorted in the same way as the embeddings
 - sample embeddings are sorted by sample name before placed in the index
 - so the ids file will be sorted in the same way.

Basic inference pipeline for evaluation of the model.
Config items:
- model: str <-  path to the model checkpoint
- gpu: boolean <- whether to use gpu
- input_files: list of str <- list of paths to input files
  NOTE: currently the positional encodings need to be offset
        back to 0 to get relative positions.

- outdir: str <- path to output directory
- batch_size: int <- batch size for inference
- num_workers: int <- number of workers for dataloader
# TODO add as needed
"""
configfile: "conf/eval_config.yaml"
config = SimpleNamespace(**config) # for dot access of config items

# segments = [163, 164, 165, 166, 167]
# file format for now 1KG.data.seg.N.pos_encoded where N is the segment number
segments = [x.split('.')[-2] for x in config.input_files]

samples = [line.rstrip() for line in open(config.samples_list)]

rule all:
  input:
    f"{config.outdir}/jaccard_similarities/jaccard.by_sample.pdf",
    f"{config.outdir}/jaccard_similarities/jaccard.similarities.png",


rule EncodeSequences:
  """
  Use model to encode seuqences into embedding vectors
  """
  input:
    expand("{input_file}", input_file=config.input_files)
  params:
    gpu = "--gpu" if config.gpu else ""
  output:
    f"{config.outdir}/embeddings.txt"
  threads:
    config.num_workers
  shell:
    f"""
    echo {config.outdir}
    python encode_samples.py \\
    --encoder {config.model} \\
    --output {{output}} \\
    --files {{input}} \\
    --batch-size {config.batch_size} \\
    {{params.gpu}} \\
    --num-workers {config.num_workers}
    """

rule MakeEmbeddingIndex:
  """
  Put embedding vectors into faiss IP index
  """
  input:
    rules.EncodeSequences.output
  output:
    index = f"{config.outdir}/faiss_encoded/index.{{segment}}.faiss",
    ids = f"{config.outdir}/faiss_encoded/ids.{{segment}}.txt",
  threads:
    1
  conda:
    "envs/faiss.yaml"
  shell:
    f"""
    python testing/index_from_encodings.py \\
    --segment {{wildcards.segment}} \\
    --embeddings {{input}} \\
    --outdir {config.outdir}/faiss_encoded
    """

rule MakeGoldStandardIndex:
  """
  Make a faiss index from the unencoded genotype vectors
  L2 index for now.
  """

  input:
    # TODO put this in config
    f"/data/segments/1KG.data.seg.{{segment}}.encoded"
  output:
    index = f"{config.outdir}/faiss_gold/index.{{segment}}.faiss",
    ids = f"{config.outdir}/faiss_gold/ids.{{segment}}.txt",
  threads:
    1
  conda:
    "envs/faiss.yaml"
  shell:
    f"""
    python testing/index_from_genotypes.py \\
    --file {{input}} \\
    --outdir {config.outdir}/faiss_gold
    """

rule QueryEmbeddingIndex:
  """
  Query the embedding index for the nearest neighbors
  """
  input:
    embeddings = rules.EncodeSequences.output,
    index = rules.MakeEmbeddingIndex.output.index,
    ids = rules.MakeEmbeddingIndex.output.ids
  output: 
    f"{config.outdir}/embedding_queries/queries.{{segment}}.txt"
  threads:
    1
  conda:
    "envs/faiss.yaml"
  shell:
    f"""
    python testing/query_faiss_encodings.py \\
    --encodings {config.outdir}/embeddings.txt \\
    --index {{input.index}} \\
    --ids {{input.ids}} \\
    --output {{output}} \\
    --k {config.num_neighbors}
    """

rule QueryGoldStandardIndex:
  """
  Query the gold standard index for the nearest neighbors
  """
  input:
    index = rules.MakeGoldStandardIndex.output.index,
    ids = rules.MakeGoldStandardIndex.output.ids
  params:
    query_file = '/data/segments/1KG.data.seg.{segment}.encoded'
  output: 
    f"{config.outdir}/gold_queries/queries.{{segment}}.txt"
  threads:
    1
  conda:
    "envs/faiss.yaml"
  shell:
    f"""
    python testing/query_faiss_genotypes.py \\
    --queries {{params.query_file}} \\
    --index {{input.index}} \\
    --ids {{input.ids}} \\
    --output {{output}} \\
    --k {config.num_neighbors}
    """

rule ComputeJaccardSimilarity:
  """
  Compute the Jaccard similarity between the embedding query results and the gold standard query results
  """
  input:
    embedding_queries = rules.QueryEmbeddingIndex.output,
    gold_queries = rules.QueryGoldStandardIndex.output
  output:
    f"{config.outdir}/jaccard_similarities/jaccard.similarity.{{segment}}.txt"
  threads:
    1
  shell:
    f"""
    python testing/compute_jaccard_similarity.py \\
    --embedding-queries {{input.embedding_queries}} \\
    --gold-queries {{input.gold_queries}} \\
    --output {{output}}
    """

rule VisualizeJaccardSimilarity:
  """
  Visualize the Jaccard similarity results
  """
  input:
    expand(rules.ComputeJaccardSimilarity.output, segment=segments)
  output:
    f"{config.outdir}/jaccard_similarities/jaccard.similarities.png",
  threads:
    1
  conda: 
    "envs/matplotlib.yaml"
  shell:
    f"""
    python testing/visualize_jaccard_similarity.py \\
    --input {{input}} \\
    --output {{output}}
    """

rule VisualizeJaccardBySample:
  """
  Sample level view of the Jaccard similarity results
  """
  input:
    expand(rules.ComputeJaccardSimilarity.output, segment=segments)
  output:
    f"{config.outdir}/jaccard_similarities/jaccard.by_sample.pdf",
  threads:
    1
  conda:
    "envs/matplotlib.yaml"
  shell:
    f"""
    python testing/visualize_jaccard_by_sample.py \\
    --input {{input}} \\
    --output {{output}}
    """
  
