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


def jaccard_similarity(set1: set, set2: set):
    return len(set1.intersection(set2)) / len(set1.union(set2))


def main(embedding_queries, gold_queries, output):
    """
    Query format is:
        <sample_id>\t<query_sample_id>:<score>\t ...
    """

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

    similarities = []
    assert gold_samples == embedding_samples
    for i in range(len(gold_samples)):
        gold_sample = gold_samples[i]
        embedding_sample = embedding_samples[i]
        assert gold_sample == embedding_sample

        similarities.append(
            jaccard_similarity(gold_results[i], embedding_results[i]),
        )

    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as f:
        for sample_id, similarity in zip(gold_samples, similarities):
            f.write(f"{sample_id}\t{similarity}\n")


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
