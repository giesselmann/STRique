# Tests

After installation, the successful setup can be tested with sample data included in the STRique repository. Test the pipeline with the following commands in the cloned repository:

```
cd ~/src/STRique
python3 scripts/STRique_test.py
python3 scripts/STRique.py index data/ > data/reads.fofn
cat data/c9orf72.sam | python3 scripts/STRique.py count ./data/reads.fofn ./models/r9_4_450bps.model ./configs/repeat_config.tsv --config ./configs/STRique.json
```

You should see output similar to

```
ID      target strand count score_prefix score_suffix log_p offset ticks mod
ce47b364-ed6e-4409-808a-1041c0b5aac2 c9orf72 - 735 6.3155927807600545 6.031860427335506 -119860.52066647023 1633 40758 -
```

Run time <1 min for this dataset on a typical desktop computer.
