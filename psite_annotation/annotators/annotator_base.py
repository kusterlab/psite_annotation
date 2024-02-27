from typing import Protocol
from functools import wraps

import pandas as pd


class Annotator(Protocol):
    """Protocol for annotating a pandas dataframe."""

    def load_annotations(self) -> None:
        """Load annotations, e.g. from a csv file."""

    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add annotation columns to dataframe.

        Args:
            df: pandas dataframe to be annotated

        """


class MissingColumnsError(Exception):
    pass


def check_columns(columns):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Assuming the DataFrame is the first argument
            df = args[1] if args and len(args) > 0 else kwargs.get("df")

            if not isinstance(df, pd.DataFrame):
                raise ValueError("The first argument should be a pandas DataFrame.")

            missing_columns = [col for col in columns if col not in df.columns]

            if missing_columns:
                instructions = ""
                if "Site sequence context" in missing_columns:
                    instructions = "You can add this column by running addSiteSequenceContext()."
                elif len(set(["Matched proteins", "Start positions", "End positions", "Site positions"]).intersection(missing_columns)) > 0:
                    instructions = "You can add these columns by running addPeptideAndPsitePositions()."
                raise MissingColumnsError(
                    f"Missing columns: {', '.join(missing_columns)}. {instructions}"
                )

            return func(*args, **kwargs)

        wrapper.__doc__ = f"""{func.__doc__}

        Required columns:
            {', '.join(map(lambda x: f':code:`{x}`', columns))}"""
        
        return wrapper

    return decorator
