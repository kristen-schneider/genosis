import torch
from transformers import LongformerTokenizerFast
from pprint import pprint

# dummy text with:
# - variable length input
# - cls, sep, pad
text =[
        "[CLS]The quick brown [MASK] jumps over the lazy dog.[SEP]Hello world!",
        "[CLS]Then he [MASK] the [MASK] and [MASK] it.[SEP]Yep that's it!",
        ]

# load tokenizer
# pretrained for longformer tokenizer
tokenizer = LongformerTokenizerFast.from_pretrained("allenai/longformer-base-4096")
