# Deep Learning for Image and Text: Classification and Generative Models

Collection of deep learning experiments covering **image classification**, **text classification**, and **image-to-image translation**, implemented in PyTorch. The repository benchmarks architectures ranging from convolutional networks to recurrent models and pretrained transformers, and extends to a conditional GAN for image dehazing.

> **Status:** Coursework and research notebooks — CIMAT, M.Sc. in Computer Science.

---

## Overview

This repository is organized into two research tracks:

- **Classification** — supervised models applied to image and text tasks:
  - A **CNN** trained from scratch on the Kaggle *Cats vs Dogs* dataset.
  - A family of **RNN / GRU / LSTM** models for emotion classification from short English messages.
  - A **fine-tuned GPT-2 (124M)** classifier over the same emotion dataset, benchmarked against the recurrent baseline.

- **Generative** — conditional image-to-image translation:
  - A **Pix2Pix** framework with a **U-Net VGG16** generator for single-image **dehazing** on the RESIDE-6K dataset.

The goal is to compare architectural families under a consistent protocol and to explore how large pretrained models transfer to constrained classification tasks and to dense generative tasks.

---

## Repository Structure

```
.
├── classification
│   ├── cnn
│   │   └── cnn_classification.ipynb
│   ├── rnn
│   │   └── rnn_gru_lstm_classification.ipynb
│   └── gpt2
│       └── gpt2_finetuning_classification.ipynb
├── generative
│   └── pix2pix_image_translation.ipynb
├── requirements.txt
└── README.md
```

---

## Data Sources

- **Cats vs Dogs** — Kaggle image dataset (~25,000 images) used by the CNN notebook.
- **Emotion Dataset for Emotion Recognition Tasks** — Kaggle (`parulpandey/emotion-dataset`), used by the RNN and GPT-2 notebooks. Six classes: *sadness, joy, love, anger, fear, surprise*. Predefined `train.csv`, `validation.csv`, `test.csv` splits are used as-is.
- **RESIDE-6K + BSD500** — Clean outdoor images plus synthetic haze pairs, used by the Pix2Pix dehazing notebook.

---

## Pipeline Overview

### 1. Image Classification — CNN

- Dataset loaded via `torchvision.datasets.ImageFolder`.
- Preprocessing: resize to `224 × 224`, tensor conversion, and normalization to `[-1, 1]` over the three RGB channels.
- Custom CNN parameterized by input image size and number of output classes.
- Trained with cross-entropy loss and SGD.

### 2. Text Classification — RNN / GRU / LSTM

- Tokenization with `tiktoken` using the GPT-2 encoding (`vocab_size ≈ 50,257`), fixed sequence length of 64 tokens with zero padding.
- PyTorch `Dataset` wrapping the encoded emotion corpus.
- **Part I:** BiLSTM baseline — a bidirectional LSTM with embedding, 2 recurrent layers, dropout, and a linear head over 6 emotions. Includes confusion matrix and per-class metrics.
- **Part II:** Architecture comparison — SimpleRNN (tanh), GRU, and LSTM, all bidirectional, trained 10 times each with seeds `1000–1009` on fixed train/val/test splits. Metrics (accuracy, precision, recall, F1 macro) are reported as mean ± std with box plots.
- **Part III:** Statistical significance — paired t-test (`scipy.stats.ttest_rel`) on F1 macro across the 10 runs to determine whether performance differences are significant at α = 0.05.

### 3. Text Classification — GPT-2 Fine-tuning

- Backbone: **GPT-2 small (124M)** loaded following the setup from Sebastian Raschka's *Build a Large Language Model From Scratch* (Chapter 6).
- Same tokenizer as the RNN notebook (`tiktoken` GPT-2) with `<|endoftext|>` used as the padding token.
- Classification head: `nn.Linear(768, 6)` applied to the last non-padding token's hidden state.
- Fine-tuning restricted to the last transformer block, the final `LayerNorm`, and the classification head; the rest of the pretrained weights are frozen.
- Direct benchmark against the recurrent models on the same train / validation / test partitions.

### 4. Image-to-Image Translation — Pix2Pix Dehazing

- Task: single-image dehazing formulated as paired image-to-image translation.
- Generator: **U-Net with a VGG16 encoder** (pretrained ImageNet weights).
- Discriminator: PatchGAN as in the original Pix2Pix formulation.
- Loss: adversarial (BCE) + L1 reconstruction term.
- Data augmentation with `albumentations`; training on Google Colab (T4 GPU), checkpoints persisted to Google Drive.
- Synthetic haze is generated on-the-fly at multiple densities from clean BSD500 / RESIDE-6K images.

---

## Evaluation

Each notebook reports its own metrics:

- **CNN:** training / validation accuracy and loss curves on Cats vs Dogs.
- **RNN / GRU / LSTM:** accuracy, precision, recall, and F1 macro averaged over 10 independent runs, with box plots and paired t-tests for statistical significance.
- **GPT-2:** accuracy and per-class F1 on the emotion test set, benchmarked against the recurrent baselines on the same splits.
- **Pix2Pix:** qualitative dehazing samples and quantitative PSNR / SSIM on held-out validation pairs.

---

## Dependencies

- torch, torchvision
- numpy, pandas, scipy, scikit-learn
- matplotlib, seaborn
- tiktoken, kagglehub
- albumentations, Pillow
- tqdm

---

## Execution Pipeline

1. Install dependencies from `requirements.txt`.
2. For the emotion notebooks, configure Kaggle credentials so that `kagglehub` can download `parulpandey/emotion-dataset`.
3. For the CNN notebook, place the Cats vs Dogs dataset under the paths expected by `ImageFolder`.
4. For the Pix2Pix notebook, mount Google Drive and provide `kaggle.json` to allow RESIDE-6K download; checkpoints are saved to `MyDrive/MLTarea3/checkpoints`.
5. Run each notebook independently — they are self-contained.

---

## Key Contributions

- Unified PyTorch implementation of CNN, RNN, GRU, LSTM, and GPT-2 classifiers under comparable training protocols.
- Rigorous multi-seed (×10) evaluation of recurrent architectures with paired t-tests for statistical significance.
- Head-to-head comparison between recurrent baselines and a fine-tuned pretrained transformer on the same emotion classification benchmark.
- Pix2Pix formulation of image dehazing with a VGG16-based U-Net generator and PatchGAN discriminator.
- Reproducible pipelines with fixed seeds, deterministic data splits, and Drive-persisted checkpoints for long-running experiments.

---

## Notes

- The emotion dataset splits are never mixed across experiments to keep RNN and GPT-2 results directly comparable.
- GPT-2 fine-tuning follows Raschka Ch. 6 conventions (last-token pooling, partial unfreezing).
- Pix2Pix training assumes a single T4 GPU; batch size and image resolution are tuned accordingly.

---

## Author

**Kaleb Aguilar**

Physics Engineer
M.Sc. Student in Computer Science at CIMAT

LinkedIn:
https://linkedin.com/in/kaleb-aguilar-70971129a
