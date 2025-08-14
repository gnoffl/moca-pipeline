# ==================================================================
# File: moca/steps/step10_browser_viz_specialized.py (Consolidated Final Version)
# ==================================================================
import os
import pandas as pd
import json
import re
import numpy as np
import base64
from pathlib import Path

# --- Helper Functions ---

def _load_gff(filepath):
    """Loads a GFF file, parsing attributes to link features for detailed models."""
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        df['id'] = df['attributes'].str.extract(r'ID=([^;]+)', expand=False)
        df['parent_id'] = df['attributes'].str.extract(r'Parent=([^;]+)', expand=False)
        return df
    except Exception as e:
        print(f"Warning: Could not load or parse GFF file {filepath}: {e}")
        return None

def image_to_base64(path):
    """Converts an image file to a base64 string for HTML embedding."""
    if not path or not Path(path).exists():
        return "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs="
    try:
        with open(path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode()
        return f"data:image/png;base64,{encoded_string}"
    except Exception as e:
        print(f"Warning: Could not encode image {path} to base64. {e}")
        return "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs="

def get_top_match(motif_id, comparison_dir, db_name):
    """Finds the top comparison hit for a motif from a specific database."""
    comp_file = Path(comparison_dir) / db_name / "comparison_results.tsv"
    if not comp_file.exists():
        return "N/A"
    try:
        df = pd.read_csv(comp_file, sep='\t')
        top_hit = df[df['Query_ID'] == motif_id].sort_values('PCC', ascending=False).iloc[0]
        return f"{top_hit['DB_Name']} (PCC: {top_hit['PCC']:.2f})"
    except (IndexError, FileNotFoundError, KeyError):
        return "No match"

def generate_html_report(data, output_path):
    """Generates the final HTML report file with visualization and table."""
    json_data = json.dumps(data)
    
    db_headers = "".join([f"<th>Top {db_name} Match</th>\n" for db_name in data['db_names']])

    table_rows_html = ""
    for motif in data['motifs_table']:
        db_matches_html = "".join([f"<td>{motif['matches'].get(db_name, 'N/A')}</td>\n" for db_name in data['db_names']])
        table_rows_html += f"""
            <tr>
                <td>{motif['id']}</td>
                <td class="logo-cell"><img src="{motif['logo_fwd']}" alt="Fwd Logo"></td>
                <td class="logo-cell"><img src="{motif['logo_rev']}" alt="Rev Logo"></td>
                <td class="plot-cell"><img src="{motif['range_plot']}" alt="Range Plot"></td>
                <td>{motif['importance']:.4f}</td>
                <td>{motif['position']}</td>
                <td>{motif['distance']}</td>
                <td>{motif['strand']}</td>
                {db_matches_html}
            </tr>"""

    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>MOCA Report: {data['gene_id']}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 2em; background-color: #fdfdfd; }}
        h1, h2, h3 {{ color: #333; }}
        #chart-container {{ width: 100%; overflow-x: auto; border: 1px solid #ccc; background-color: #f9f9f9; margin-bottom: 2em;}}
        .tooltip {{ position: absolute; visibility: hidden; background: rgba(0, 0, 0, 0.8); color: #fff; padding: 8px; border-radius: 4px; font-size: 12px; pointer-events: none; z-index: 10; white-space: pre; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; font-size: 12px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; vertical-align: middle; }}
        th {{ background-color: #f2f2f2; }}
        img {{ max-width: 100%; height: auto; }}
        .logo-cell {{ width: 150px; }}
        .plot-cell {{ width: 250px; }}
        .gene-backbone {{ stroke: #555; stroke-width: 2; }}
        .feature {{ cursor: pointer; }}
        .exon {{ fill: #666666; stroke: #444444; }}
        .CDS {{ fill: #666666; stroke: #444444; }}
        .five_prime_UTR, .three_prime_UTR {{ fill: #a0a0a0; stroke: #888888; }}
        .motif {{ cursor: pointer; }}
        .axis-line, .axis-tick {{ stroke: #aaa; stroke-width: 1; }}
        .axis-label {{ font-size: 10px; fill: #555; }}
    </style>
</head>
<body>
    <h1>MOCA Analysis Report</h1>
    <h2>Gene: {data['gene_id']}</h2>
    <h4>Sequence ID: {data['chr']} | Strand: {data['gene']['strand']} | Annotated Motifs: {data['total_motifs']}</h4>
    
    <h3>Interactive Gene Model</h3>
    <div id="chart-container"></div>
    <div id="tooltip" class="tooltip"></div>

    <h3>Annotated Motif Details</h3>
    <table>
        <thead>
            <tr>
                <th>Motif ID</th>
                <th>Forward Logo</th>
                <th>Reverse Logo</th>
                <th>Positional Density</th>
                <th>Importance Score</th>
                <th>Position (in sequence)</th>
                <th>Distance to Border</th>
                <th>Match Strand</th>
                {db_headers}
            </tr>
        </thead>
        <tbody>
            {table_rows_html}
        </tbody>
    </table>

    <script>
        const data = {json_data};
        document.addEventListener("DOMContentLoaded", function() {{
            const container = document.getElementById('chart-container');
            const margin = {{ top: 40, right: 30, bottom: 40, left: 30 }};
            const width = Math.max(container.clientWidth, 1200) - margin.left - margin.right;
            const height = 120 + (data.max_track * 25);
            const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
            svg.setAttribute("width", width + margin.left + margin.right);
            svg.setAttribute("height", height + margin.top + margin.bottom);
            const g = document.createElementNS("http://www.w3.org/2000/svg", "g");
            g.setAttribute("transform", `translate(${{margin.left}}, ${{margin.top}})`);
            svg.appendChild(g);
            container.appendChild(svg);
            const scale = (pos) => ((pos - data.view_start) / (data.view_end - data.view_start)) * width;
            const tooltip = document.getElementById('tooltip');

            // Draw gene model
            const geneY = 40;
            const gene = data.gene;
            const backbone = document.createElementNS("http://www.w3.org/2000/svg", "line");
            backbone.setAttribute("x1", scale(gene.start));
            backbone.setAttribute("y1", geneY);
            backbone.setAttribute("x2", scale(gene.end));
            backbone.setAttribute("y2", geneY);
            backbone.classList.add("gene-backbone");
            g.appendChild(backbone);
            
            // Draw gene features (Exons, UTRs, etc.)
            data.features.forEach(feature => {{
                const featureRect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
                const x = scale(feature.start);
                const featureWidth = scale(feature.end) - x;
                featureRect.setAttribute("x", x);
                featureRect.setAttribute("width", Math.max(1, featureWidth));
                featureRect.setAttribute("class", `feature ${{feature.type}}`);
                const featureHeight = (feature.type === 'CDS' || feature.type === 'exon') ? 15 : 10;
                featureRect.setAttribute("y", geneY - featureHeight / 2);
                featureRect.setAttribute("height", featureHeight);
                g.appendChild(featureRect);
            }});

            // Draw labels
            data.gene_labels.forEach(label => {{
                const labelText = document.createElementNS("http://www.w3.org/2000/svg", "text");
                labelText.setAttribute("x", scale(label.pos));
                labelText.setAttribute("y", geneY - 20);
                labelText.setAttribute("text-anchor", "middle");
                labelText.setAttribute("font-size", "12px");
                labelText.setAttribute("font-weight", "bold");
                labelText.textContent = label.text;
                g.appendChild(labelText);
            }});
            
            // Draw motifs
            data.motifs_viz.forEach(motif => {{
                const motifGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
                const motifY = 65 + (motif.track * 20);
                const x1 = scale(motif.start);
                const x2 = scale(motif.end);
                const motifPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
                const d = `M${{x1}},${{motifY}} L${{x2}},${{motifY}}`;
                motifPath.setAttribute("d", d);
                motifPath.setAttribute("stroke", motif.color);
                motifPath.setAttribute("stroke-width", "5");
                motifGroup.appendChild(motifPath);
                
                const arrow = document.createElementNS("http://www.w3.org/2000/svg", "path");
                const arrowD = motif.strand === '+' ? `M${{x2-5}},${{motifY-4}} L${{x2}},${{motifY}} L${{x2-5}},${{motifY+4}}` : `M${{x1+5}},${{motifY-4}} L${{x1}},${{motifY}} L${{x1+5}},${{motifY+4}}`;
                arrow.setAttribute("d", arrowD);
                arrow.setAttribute("fill", motif.color);
                motifGroup.appendChild(arrow);
                g.appendChild(motifGroup);

                motifGroup.addEventListener('mouseover', (event) => {{
                    tooltip.style.visibility = 'visible';
                    tooltip.innerHTML = `<b>Motif:</b> ${{motif.id}}<br><b>Score:</b> ${{motif.score.toFixed(2)}}<br><b>Position:</b> ${{motif.start}} - ${{motif.end}}`;
                }});
                motifGroup.addEventListener('mousemove', (event) => {{
                    tooltip.style.top = (event.pageY + 15) + 'px';
                    tooltip.style.left = (event.pageX + 15) + 'px';
                }});
                motifGroup.addEventListener('mouseout', () => {{ tooltip.style.visibility = 'hidden'; }});
            }});

            // Draw coordinate axis
            const axisY = height - 15;
            const axisLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
            axisLine.setAttribute("x1", 0);
            axisLine.setAttribute("y1", axisY);
            axisLine.setAttribute("x2", width);
            axisLine.setAttribute("y2", axisY);
            axisLine.classList.add("axis-line");
            g.appendChild(axisLine);
            const numTicks = 10;
            for (let i = 0; i <= numTicks; i++) {{
                const pos = data.view_start + i * (data.view_end - data.view_start) / numTicks;
                const x = scale(pos);
                const tick = document.createElementNS("http://www.w3.org/2000/svg", "line");
                tick.setAttribute("x1", x);
                tick.setAttribute("y1", axisY);
                tick.setAttribute("x2", x);
                tick.setAttribute("y2", axisY + 5);
                tick.classList.add("axis-tick");
                g.appendChild(tick);
                const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
                label.setAttribute("x", x);
                label.setAttribute("y", axisY + 15);
                label.setAttribute("text-anchor", "middle");
                label.classList.add("axis-label");
                label.textContent = Math.round(pos);
                g.appendChild(label);
            }}
        }});
    </script>
</body>
</html>"""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_template)


def run(config, common_settings):
    """
    Generates a feature-rich report for each target sequence, consolidating all analysis results.
    """
    print("--- Running Full-Featured Browser Visualization ---")
    
    # --- 1. Configuration & Path Setup ---
    output_dir = config['output_dir']
    target_gene_list_file = config.get('target_gene_list')
    
    annotation_dir = common_settings.get('annotation', {}).get('output_dir')
    viz_dir = common_settings.get('visualization', {}).get('output_dir')
    ranging_dir = common_settings.get('ranging', {}).get('output_dir')
    importance_dir = common_settings.get('importance', {}).get('output_dir')
    comparison_dir = common_settings.get('comparison', {}).get('output_dir')
    ref_gff_file = common_settings.get('annotation', {}).get('reference_gff')

    filter_method = common_settings.get('annotation', {}).get('filter_method', 'minmax')
    annotated_motifs_file = Path(annotation_dir) / f"annotated_motifs_{filter_method}.csv" if annotation_dir else None
    importance_file = Path(importance_dir) / f"importance_scores.csv" if importance_dir else None

    # --- 2. Load Core Data ---
    if not annotated_motifs_file or not annotated_motifs_file.exists():
        print(f"Error: Annotated motifs file not found at '{annotated_motifs_file}'. Skipping.")
        return

    print("Loading core data...")
    annotated_df = pd.read_csv(annotated_motifs_file)
    gff_df = _load_gff(ref_gff_file) if ref_gff_file and Path(ref_gff_file).exists() else None
    importance_df = pd.read_csv(importance_file) if importance_file and importance_file.exists() else None

    with open(target_gene_list_file, 'r') as f:
        target_sequence_ids = [line.strip() for line in f if line.strip()]

    db_names = [d.name for d in Path(comparison_dir).iterdir() if d.is_dir()] if comparison_dir and Path(comparison_dir).exists() else []

    # --- 3. Main Loop: Generate a Report for Each Sequence ---
    print(f"\nFound {len(target_sequence_ids)} target sequences to visualize.")
    for seq_id in target_sequence_ids:
        print(f"--- Generating report for sequence: {seq_id} ---")

        target_motifs_df = annotated_df[annotated_df['sequence_id'] == seq_id].copy()
        if target_motifs_df.empty:
            print(f"Warning: No annotated motifs found for '{seq_id}'. Skipping report.")
            continue

        simple_gene_id = target_motifs_df['gene_id'].iloc[0]

        # --- 4. Gene Model and Feature Definition ---
        gene_info = gff_df[(gff_df['seqid'] == seq_id) & (gff_df['type'] == 'gene')] if gff_df is not None else pd.DataFrame()
        features = []
        gene_labels = []

        if gene_info.empty:
            print(f"--> No gene found in GFF. Generating a dummy model for {seq_id}.")
            gene_start, gene_end, gene_strand = 1000, 2020, '+'
            features = [{'start': 1000, 'end': 2020, 'type': 'exon'}]
            gene_labels = [{'pos': (gene_start + gene_end) / 2, 'text': 'Gene Region'}]
        else:
            gene_start, gene_end, gene_strand = int(gene_info['start'].iloc[0]), int(gene_info['end'].iloc[0]), gene_info['strand'].iloc[0]
            gene_gff_id = gene_info['id'].iloc[0]
            all_transcripts = gff_df[(gff_df['parent_id'] == gene_gff_id) & (gff_df['type'] == 'mRNA')]
            if not all_transcripts.empty:
                all_transcripts['length'] = all_transcripts['end'] - all_transcripts['start']
                longest_transcript = all_transcripts.loc[all_transcripts['length'].idxmax()]
                feature_types = ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
                feature_info = gff_df[(gff_df['parent_id'] == longest_transcript['id']) & (gff_df['type'].isin(feature_types))]
                features = [{'start': int(r['start']), 'end': int(r['end']), 'type': r['type']} for _, r in feature_info.iterrows()]

        view_start, view_end = 1, 3020

        # --- 5. Motif Layout and Color ---
        target_motifs_df = target_motifs_df.sort_values('mstart').reset_index(drop=True)
        target_motifs_df['track'] = 0
        track_ends = {0: -1}
        for i, row in target_motifs_df.iterrows():
            placed = False
            for track, end_pos in sorted(track_ends.items()):
                if row['mstart'] > end_pos + 5:
                    target_motifs_df.loc[i, 'track'] = track
                    track_ends[track] = row['mend']
                    placed = True
                    break
            if not placed:
                new_track = max(track_ends.keys()) + 1
                target_motifs_df.loc[i, 'track'] = new_track
                track_ends[new_track] = row['mend']
        
        unique_motifs = target_motifs_df['motif'].unique()
        color_map = {motif: f"hsl({(i * 360 / len(unique_motifs)) % 360}, 90%, 45%)" for i, motif in enumerate(unique_motifs)}
        target_motifs_df['color'] = target_motifs_df['motif'].map(color_map)
        
        # --- 6. Assemble Data for Report ---
        report_data = {
            'gene_id': simple_gene_id,
            'chr': seq_id,
            'view_start': view_start, 'view_end': view_end,
            'gene': {'start': gene_start, 'end': gene_end, 'strand': gene_strand},
            'features': features,
            'gene_labels': gene_labels,
            'motifs_viz': [{'id': r['motif'], 'start': r['mstart'], 'end': r['mend'], 'score': r['score'], 'track': r['track'], 'color': r['color'], 'strand': r['strand']} for _, r in target_motifs_df.iterrows()],
            'max_track': int(target_motifs_df['track'].max()) if not target_motifs_df.empty else 0,
            'total_motifs': len(target_motifs_df),
            'db_names': db_names,
            'motifs_table': []
        }
        
        for _, row in target_motifs_df.iterrows():
            motif_id = row['motif']
            is_reverse = 'R_' in motif_id
            base_motif_id = motif_id.replace('R_', 'F_')
            motif_short_name_match = re.search(r'(p\d+m\d+)', base_motif_id)
            motif_short_name = motif_short_name_match.group(1) if motif_short_name_match else ""

            logo_fwd_path = Path(viz_dir) / f"{base_motif_id}.png" if viz_dir else None
            logo_rev_path = Path(viz_dir) / f"{base_motif_id.replace('F_', 'R_')}.png" if viz_dir else None
            
            range_plot_path = None
            if ranging_dir and motif_short_name:
                species_tag = common_settings.get('species_tag', 'unk')
                model_tag = common_settings.get('model_tag', 'm0')
                plot_filename = f"epm_{species_tag}_{model_tag}_{motif_short_name}_density.png"
                range_plot_path = Path(ranging_dir) / 'distribution_plots' / plot_filename
            
            importance_score = 0.0
            if importance_df is not None:
                score_row = importance_df[importance_df['motif'] == motif_id]
                if not score_row.empty:
                    if 'importance_score' in score_row.columns:
                         importance_score = score_row['importance_score'].iloc[0]
                    elif 'contrib_score_sum' in score_row.columns:
                         importance_score = score_row['contrib_score_sum'].iloc[0]

            matches = {db_name: get_top_match(motif_id, comparison_dir, db_name) for db_name in db_names}

            report_data['motifs_table'].append({
                'id': motif_id,
                'logo_fwd': image_to_base64(logo_rev_path) if is_reverse else image_to_base64(logo_fwd_path),
                'logo_rev': image_to_base64(logo_fwd_path) if is_reverse else image_to_base64(logo_rev_path),
                'range_plot': image_to_base64(range_plot_path),
                'importance': importance_score,
                'position': f"{row['mstart']}-{row['mend']}",
                'distance': f"{row['dist_transc_border']}bp ({row['region']})",
                'strand': row['strand'],
                'matches': matches
            })

        # --- 7. Generate and Save HTML with a Unique Filename ---
        unique_suffix = ""
        # Find a unique part like 'mutations_018' from the long sequence ID
        match = re.search(r'(mutations_\d+)', seq_id)
        if match:
            unique_suffix = match.group(1)
        else:
            # Fallback if pattern isn't found: use the last 20 sanitized characters of the unique seq_id
            sanitized_end = re.sub(r'[^a-zA-Z0-9_-]', '_', seq_id[-20:])
            unique_suffix = sanitized_end

        # Combine the base gene ID with the unique suffix
        final_filename_base = f"{simple_gene_id}_{unique_suffix}"
        output_path = Path(output_dir) / f"{final_filename_base}_moca_report.html"
        generate_html_report(report_data, output_path)
        print(f"--> Interactive HTML report saved to: {output_path}")

    print("\nBrowser visualization step successfully completed.")
