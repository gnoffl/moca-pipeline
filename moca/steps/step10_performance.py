import os
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from tqdm import tqdm

def _clean_gene_id(series):
    """Standardizes gene IDs for merging by lowercasing and removing suffixes."""
    return series.str.lower().str.replace(r"[-']+", "", regex=True).str.replace(r"\..*", "", regex=True)

def run(config, common_settings):
    """
    Main function for the performance evaluation step.
    Evaluates motif predictiveness and performs functional enrichment analysis.
    """
    # --- 1. Configuration ---
    output_dir = config['output_dir']
    
    # Input files
    annotation_dir = common_settings.get('annotation', {}).get('output_dir')
    annotated_motifs_file = os.path.join(annotation_dir, "annotated_motifs_q1q9.csv")
    
    predictions_file = config.get('prediction_probabilities_csv')
    expression_file = config.get('true_expression_classes_csv')
    functional_annot_file = config.get('functional_annotation_txt')

    if not all(os.path.exists(f) for f in [annotated_motifs_file, expression_file, functional_annot_file]):
        print("Error: One or more required input files not found. Skipping performance analysis.")
        return

    # --- 2. Load and Prepare Data ---
    print("Loading data for performance analysis...")
    annotated_motifs = pd.read_csv(annotated_motifs_file, low_memory=False)
    expression_df = pd.read_csv(expression_file)
    func_annot_df = pd.read_csv(functional_annot_file, sep='\t', quotechar='"')

    # Handle optional prediction probabilities file
    if predictions_file and os.path.exists(predictions_file):
        preds_df = pd.read_csv(predictions_file)
        preds_df.columns = ['gene_id', 'prob']
    else:
        print("Warning: Prediction probabilities file not found. Generating dummy data.")
        unique_genes = annotated_motifs['gene_id'].unique()
        preds_df = pd.DataFrame({
            'gene_id': unique_genes,
            'prob': np.random.randint(0, 2, size=len(unique_genes))
        })

    # --- 3. Clean IDs and Merge DataFrames ---
    print("Standardizing gene IDs and merging data sources...")
    annotated_motifs['gene_id_clean'] = _clean_gene_id(annotated_motifs['gene_id'])
    preds_df['gene_id_clean'] = _clean_gene_id(preds_df['gene_id'])
    expression_df['gene_id_clean'] = _clean_gene_id(expression_df.columns[0]) # Assuming gene ID is the first column
    func_annot_df['gene_id_clean'] = _clean_gene_id(func_annot_df['IDENTIFIER'])

    # Rename columns for clarity
    expression_df.rename(columns={expression_df.columns[1]: 'class'}, inplace=True)
    
    # Merge all data
    merged = pd.merge(annotated_motifs, preds_df, on='gene_id_clean', how='inner')
    merged = pd.merge(merged, expression_df[['gene_id_clean', 'class']], on='gene_id_clean', how='inner')
    
    # --- 4. Calculate Performance Metrics ---
    print("Calculating performance metrics...")
    merged['prob_class'] = np.where(merged['prob'] > 0.5, 'high', 'low')
    merged['class'] = np.where(merged['class'] == 1, 'high', 'low')
    merged['pred_perf'] = merged['prob_class'] == merged['class']

    is_p0 = merged['epm'].str.contains('_p0m')
    is_p1 = merged['epm'].str.contains('_p1m')
    
    cond_true_high = is_p0 & (merged['class'] == 'high') & merged['pred_perf']
    cond_true_low = is_p1 & (merged['class'] == 'low') & merged['pred_perf']
    
    merged['epm_pred_perf'] = np.select(
        [cond_true_high, cond_true_low], ['TRUE_high', 'TRUE_low'], default='NA'
    )

    # --- 5. Create Contingency Tables and Statistics ---
    contingency_expr = pd.crosstab(merged['epm'], merged['class'])
    contingency_prob = pd.crosstab(merged['epm'], merged['prob_class'])
    
    perf_subset = merged[merged['epm_pred_perf'] != 'NA']
    contingency_perf = pd.crosstab(perf_subset['epm'], perf_subset['epm_pred_perf'])

    summary_df = pd.concat([contingency_expr, contingency_prob, contingency_perf], axis=1).fillna(0)
    summary_df.columns = ['expr_high', 'expr_low', 'prob_high', 'prob_low', 'perf_TRUE_high', 'perf_TRUE_low']

    summary_df['epm_TPpred'] = summary_df['perf_TRUE_high'] + summary_df['perf_TRUE_low']
    summary_df['epm_total'] = summary_df['expr_high'] + summary_df['expr_low']
    summary_df['TPR_TF'] = summary_df['epm_TPpred'] / summary_df['epm_total']

    # Save performance results
    perf_path = os.path.join(output_dir, "motif_performance_summary.csv")
    summary_df.to_csv(perf_path, float_format='%.4f')
    print(f"Motif performance summary saved to {perf_path}")

    # --- 6. Functional Enrichment Analysis ---
    print("Performing functional enrichment analysis...")
    merged_func = pd.merge(merged, func_annot_df[['gene_id_clean', 'BINCODE', 'NAME']], on='gene_id_clean', how='inner')
    
    # Create contingency table for motifs vs. functional bins
    func_contingency = pd.crosstab(merged_func['epm'], merged_func['BINCODE'])
    
    # Perform Chi-squared test
    chi2, p, dof, expected = chi2_contingency(func_contingency)
    
    # Find the most significant associations
    chi2_residuals = (func_contingency - expected) / np.sqrt(expected)
    
    # Unstack to get a series of (epm, BINCODE) pairs and their residual values
    enrichment_scores = chi2_residuals.stack().reset_index()
    enrichment_scores.columns = ['epm', 'BINCODE', 'chi2_residual']
    
    # Get the functional category name
    bin_names = func_annot_df[['BINCODE', 'NAME']].drop_duplicates()
    enrichment_scores = pd.merge(enrichment_scores, bin_names, on='BINCODE', how='left')
    
    # Sort by the strength of association and save the top results
    top_enrichment = enrichment_scores.sort_values(by='chi2_residual', ascending=False).head(100)
    
    enrich_path = os.path.join(output_dir, "top100_functional_enrichment.csv")
    top_enrichment.to_csv(enrich_path, index=False, float_format='%.3f')
    print(f"Top 100 functional enrichment results saved to {enrich_path}")

    print("Performance analysis step successfully completed.")
