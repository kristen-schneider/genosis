import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--embedding-queries",
        required=True,
        help="Path to the embedding query results",
    )
    parser.add_argument(
        "--gold-queries",
        required=True,
        help="Path to the gold standard query results",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path of output file",
    )
    return parser.parse_args()


def parse_query_line(line):
    A = line.strip().split("\t")
    sample_id = A[0]
    query_results = set()
    for query_result in A[1:]:
        query_sample_id, score = query_result.split(":")
        query_results.add(query_sample_id)
    return sample_id, query_results


def precision_recall(gold_standard: set, predictions: set):
    TP = len(gold_standard.intersection(predictions))
    FP = len(predictions.difference(gold_standard))
    FN = len(gold_standard.difference(predictions))
    print(f"{TP=}, {FP=}, {FN=}")
    return TP / (TP + FP), TP / (TP + FN)


def main(embedding_queries, gold_queries, output):
    """
    Query format is:
        <sample_id>\t<query_sample_id>:<similarity>\t ...
    """

    # load gold standard
    gold_samples = []
    gold_results = []
    with open(gold_queries) as f:
        for line in f:
            s, r = parse_query_line(line)
            gold_samples.append(s)
            gold_results.append(r)
    gold_samples, gold_results = zip(
        *sorted(zip(gold_samples, gold_results)),
    )

    # load embedding results
    embedding_samples = []
    embedding_results = []
    with open(embedding_queries) as f:
        for line in f:
            s, r = parse_query_line(line)
            embedding_samples.append(s)
            embedding_results.append(r)
    embedding_samples, embedding_results = zip(
        *sorted(zip(embedding_samples, embedding_results)),
    )

    # get precision and recall of each query
    scores = []
    assert gold_samples == embedding_samples
    for i in range(len(gold_samples)):
        gold_sample = gold_samples[i]
        embedding_sample = embedding_samples[i]
        assert gold_sample == embedding_sample

        scores.append(
            precision_recall(gold_results[i], embedding_results[i]),
        )

    # write results with format <sample_id>\t<precision>\t<recall>
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as f:
        for sample_id, (p, r) in zip(gold_samples, scores):
            f.write(f"{sample_id}\t{p}\t{r}\n")


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
