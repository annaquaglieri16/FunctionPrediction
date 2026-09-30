"""Frozen PLM embedding extraction.

Deliberately lazy-imports torch/transformers so the rest of this package
(digestion, fragment generation, tests) works without the heavy ML
dependencies installed. Install requirements-ml.txt before using this
module.
"""

from __future__ import annotations

import numpy as np

_MODEL_CACHE: dict[str, tuple] = {}

DEFAULT_MODEL = "facebook/esm2_t12_35M_UR50D"  # small ESM2 checkpoint, fast to iterate on


def _load_model(model_name: str):
    if model_name in _MODEL_CACHE:
        return _MODEL_CACHE[model_name]

    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "embed.py needs torch + transformers. Install requirements-ml.txt first."
        ) from exc

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    _MODEL_CACHE[model_name] = (tokenizer, model, torch)
    return _MODEL_CACHE[model_name]


POOLING_MODES = ("mean", "mean_with_special", "cls")


def embed_sequences(
    sequences: list[str],
    model_name: str = DEFAULT_MODEL,
    pooling: str = "mean",
    batch_size: int = 8,
    layer: int = -1,
) -> np.ndarray:
    """Return one fixed-length frozen embedding per input sequence.

    pooling:
      "mean" (default) -- average token embeddings over the *residue*
          positions only, excluding padding **and** the special tokens
          (`<cls>`/`<eos>`) that the ESM tokenizer adds. This matters a lot
          for short sequences: for a 6-mer the two special tokens are 2 of 8
          positions, i.e. 25% of the pooled vector would otherwise be two
          constant vectors shared by every sequence. That constant share
          scales as 1/length, which makes the pooled embedding encode length
          and turns cosine similarity into a length-similarity detector
          (measured on the notebook-06 pilot: corr(cosine, |dlength|) = -0.33,
          and a same-function AUROC pushed *below* chance).
      "mean_with_special" -- the old behaviour, kept so pre-fix results stay
          reproducible. Do not use for new runs.
      "cls" -- the first token's embedding. Note ESM-2 has no sentence-level
          pretraining objective, so `<cls>` is not trained to summarise the
          sequence; this is usually worse than masked mean pooling.

    layer: which hidden layer to pool. -1 (default) uses the final layer.
        Any other value indexes `output_hidden_states` (0 = embedding layer),
        since intermediate layers are often better for similarity tasks than
        the final layer, which is specialised for the masked-token head.
    """
    if pooling not in POOLING_MODES:
        raise ValueError(f"pooling must be one of {POOLING_MODES}")

    tokenizer, model, torch = _load_model(model_name)

    all_vecs = []
    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i + batch_size]
            enc = tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                truncation=True,
                return_special_tokens_mask=True,
            )
            special = enc.pop("special_tokens_mask")
            out = model(**enc, output_hidden_states=(layer != -1))
            hidden = out.last_hidden_state if layer == -1 else out.hidden_states[layer]

            if pooling == "cls":
                vecs = hidden[:, 0, :]
            else:
                mask = enc["attention_mask"]
                if pooling == "mean":
                    # real residues only: attended AND not a special token
                    mask = mask * (1 - special)
                mask = mask.unsqueeze(-1).to(hidden.dtype)
                summed = (hidden * mask).sum(dim=1)
                counts = mask.sum(dim=1).clamp(min=1)
                vecs = summed / counts

            all_vecs.append(vecs.cpu().numpy())

    return np.concatenate(all_vecs, axis=0)
