import sys
sys.path.insert(0, "..")
import merge_polarity

def test_merge():
    merge_polarity.merge("positive.graphml", "negative.graphml", "merged.graphml", output_summary_table="summary.tsv")