import pandas as pd

from sepsis_deescalation.pcornet_io import combine_date_time


def test_combine_date_time_text_and_integer_hhmm():
    dates = pd.Series(["2020-01-01", "2020-01-02"])
    text = pd.Series(["13:45", "00:05"])
    out = combine_date_time(dates, text)
    assert out.iloc[0] == pd.Timestamp("2020-01-01 13:45")
    assert out.iloc[1] == pd.Timestamp("2020-01-02 00:05")

    ints = pd.Series([1345, 5])
    out2 = combine_date_time(dates, ints)
    assert out2.iloc[0] == pd.Timestamp("2020-01-01 13:45")
    assert out2.iloc[1] == pd.Timestamp("2020-01-02 00:05")
