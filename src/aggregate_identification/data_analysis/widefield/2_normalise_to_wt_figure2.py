"""Normalize aggregate-positive percentages to WT for Figure 2 widefield data.

This script is a Figure 2-specific follow-on step for
`1_aggregate_count_plate-map.py`. It reads the exported `count_melt` table,
filters to percentage-positive values, normalizes each observation to the WT
reference condition, and saves the normalized dataframe for downstream plotting.
"""

import os
import pandas as pd


def add_treatment_components(df: pd.DataFrame) -> pd.DataFrame:
    """Split the treatment string into component columns.

    Args:
        df: Dataframe containing a `treat` column with underscore-delimited
            treatment labels.

    Returns:
        Dataframe with `t1`, `t2`, and `t3` columns added when possible.
    """
    if 't3' in df.columns:
        return df

    split_cols = df['treat'].str.split('_', expand=True)
    if split_cols.shape[1] < 3:
        raise ValueError(
            "Expected `treat` values with at least three underscore-delimited parts to derive `t3`."
        )

    df = df.copy()
    df[['t1', 't2', 't3']] = split_cols.iloc[:, :3]
    return df


def normalise_to_wt(df: pd.DataFrame, wt_label: str = 'UD-WT') -> pd.DataFrame:
    """Normalize percentage-positive values to WT within each timepoint.

    Args:
        df: Long-format dataframe containing at least `tp`, `t3`, and `value`.
        wt_label: Value in `t3` that identifies the WT reference condition.

    Returns:
        Dataframe with an added `normed_to_WT` column.
    """
    df = df.copy()

    wt_by_tp = (
        df[df['t3'] == wt_label]
        .groupby('tp', as_index=True)['value']
        .mean()
    )

    if wt_by_tp.empty:
        raise ValueError(f"No WT rows found for label `{wt_label}`.")

    missing_tp = sorted(set(df['tp']) - set(wt_by_tp.index))
    if missing_tp:
        raise ValueError(
            f"Missing WT rows for timepoint(s): {missing_tp}. Cannot normalize these groups."
        )

    if (wt_by_tp == 0).any():
        bad_tp = wt_by_tp[wt_by_tp == 0].index.tolist()
        raise ValueError(f"WT mean is zero for timepoint(s): {bad_tp}. Division by zero is not allowed.")

    df['normed_to_WT'] = df.apply(lambda row: row['value'] / wt_by_tp.loc[row['tp']], axis=1)
    return df


if __name__ == '__main__':
    exp_num = 'Expnum-replicate'
    wt_label = 'UD-WT'

    input_folder = f'output/Figure_2/data_analysis/1_aggregate_count/{exp_num}/'
    output_folder = f'output/Figure_2/data_analysis/1_aggregate_count/{exp_num}/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    input_path = f"{input_folder}{exp_num}_count_melt.csv"
    count_melt = pd.read_csv(input_path)
    count_melt = count_melt.drop(columns=[c for c in count_melt.columns if 'Unnamed:' in c], errors='ignore')

    # Keep only percentage-positive measurements for WT normalization.
    perc_pos = count_melt[count_melt['variable'] == 'perc_pos'].copy()
    perc_pos = add_treatment_components(perc_pos)

    normed = normalise_to_wt(perc_pos, wt_label=wt_label)
    normed.to_csv(f'{output_folder}normed_amount_agg.csv', index=False)

    print(f'Saved normalized dataframe to {output_folder}normed_amount_agg.csv')
