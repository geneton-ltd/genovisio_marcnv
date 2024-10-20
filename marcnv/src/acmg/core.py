import os

SRC_DIR = os.path.abspath(os.path.dirname(__file__))

ACMG_GAIN_TSV_FILEPATH = os.path.join(SRC_DIR, "data", "acmg_gain.tsv")
ACMG_LOSS_TSV_FILEPATH = os.path.join(SRC_DIR, "data", "acmg_loss.tsv")

# ACMG criteria constants
DUPLICATION_GENES_THRESHOLDS = (35, 50)
DELETION_GENES_THRESHOLDS = (25, 35)

HI_TS_SCORES = [3]
MIN_FREQUENCY_BENIGN = 0.005
