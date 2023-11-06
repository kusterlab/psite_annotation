"""Functions to combine separator-separated strings, e.g. with commas or semicolons."""

from itertools import chain
from typing import Dict, List, Union

import pandas as pd


def combine_separated_strings(annotations: List[str], sep=";") -> str:
    """Combines a list of separator-separated strings and returns the unique, sorted, and concatenated values.

    Args:
        annotations: A list of strings that contain separator-separated values.
        sep: The separator used in the strings. Default: ";".

    Returns:
        str: The unique, sorted, and concatenated values.
    """
    annotations = _split_and_flatten(annotations, sep=sep)
    return keep_unique_and_join(annotations, sep=sep)


def _split_and_flatten(list_of_separated_strings: List[str], sep=";"):
    """Splits a list of strings that contain separator-separated values, and flattens the list.

    Args:
        list_of_separated_strings: A list of strings.
        sep: The separator used to concatenate the strings. Default: ";".

    Returns:
        str: The concatenated unique strings.

    """
    list_of_lists = map(lambda x: x.split(sep), list_of_separated_strings)

    # flatten the list of lists into a single list
    return chain.from_iterable(list_of_lists)


def keep_unique_and_join(list_of_strings: List[str], sep: str = ";") -> str:
    """Keeps unique elements of a list and join them into a single string.

    Args:
        list_of_strings: A list of strings.
        sep: The separator used to join the strings in the list. Default is ';'.

    Returns:
        str: A single string containing all unique and sorted elements from the input, separated by the given separator.
          Note that empty strings are removed.
    """
    # Remove empty strings
    list_of_strings = filter(None, list_of_strings)
    return sep.join(sorted(set(list_of_strings)))


def merge_on_separated_string(
    df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    index_col: str,
    agg_func: Union[str, Dict],
    sep: str = ";",
):
    """Merges a dataframe with a column with semicolon separated IDs with a dataframe with a single ID column.

    Args:
        df: The original dataframe to merge with the annotation dataframe.
        annotation_df: The dataframe containing annotations for the index_col in the original dataframe.
        index_col: The column name of the index column in both dataframes that contains the ID's, note that this has
            to be an actual column and not (part of) the index of the dataframe.
        agg_func: The aggregation function to apply to the merged dataframe.
        sep: The separator used in the index_col of the original dataframe to separate the IDs. Default is ";".

    Returns:
        pd.DataFrame: The merged dataframe. Note that the index of the input dataframe will be lost.

    """
    # set aggregation for the original dataframe to simply take the original value
    agg_func_df = {k: "first" for k in df.columns}
    agg_func_merged = {**agg_func_df, **agg_func}

    # discard the original index since there might be duplicated indices, e.g. due to pd.concat
    df_copy = df.reset_index(drop=True)

    # use the native 0-based index to keep track of which rows belong together after explode
    df_copy = df_copy.reset_index(drop=False)

    df_exploded = explode_on_separated_string(df_copy, index_col, sep)
    df_merged = df_exploded.merge(annotation_df, on=index_col, how="left")

    # transform integer annotations to strings so we can concatenate them after grouping.
    # first transform back to integers, because the NaN values caused in the merge by missing indices result in the
    # column to have floats as the dtype.
    integer_annotation_columns = annotation_df.select_dtypes(include=["int"]).columns
    df_merged[integer_annotation_columns] = (
        df_merged[integer_annotation_columns].fillna(0).astype(int).astype(str)
    )

    string_annotation_columns = annotation_df.select_dtypes(include=["object"]).columns
    df_merged[string_annotation_columns] = df_merged[string_annotation_columns].fillna(
        ""
    )

    # merge rows from the same row before explosion into a single row
    df_merged = df_merged.groupby("index")[df_merged.columns].agg(agg_func_merged)

    # drop the 0-based 'index' column
    df_merged = df_merged.reset_index(drop=True)

    return df_merged


def explode_on_separated_string(df: pd.DataFrame, index_col: str, sep: str = ";"):
    """Explodes a dataframe's column with semicolon separated strings into multiple rows, one for each string.

    Args:
        df: The dataframe to explode.
        index_col: The column name of the column containing the semicolon separated strings.
        sep: The separator used in the index_col to separate the strings. Default is ";".

    Returns:
        pd.DataFrame: The exploded dataframe.
    """
    df[index_col] = df[index_col].str.split(sep)
    return df.explode(index_col)
