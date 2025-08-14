import os
import pandas as pd
import json
import re
import numpy as np
import base64
from pathlib import Path

def _normalize_gene_id(series):
    """
    Removes common prefixes and transcript versions from gene IDs for consistent matching.
    e.g., 'gene:AT1G01010.1' becomes 'AT1G01010'.
    """
    if series is None:
        return None
    return series.str.replace('^gene:', '', regex=True).str.replace(r'\.\d+$', '', regex=True)

def _load_gff(filepath):
    """Loads a GFF/GTF file, parsing it into a pandas DataFrame."""
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        
        # Extract ID, Parent, and a potential gene_id from the attributes for linking
        df['id'] = df['attributes'].str.extract(r'ID=([^;]+)', expand=False)
        df['parent_id'] = df['attributes'].str.extract(r'Parent=([^;]+)', expand=False)
        df['gene_id_attr'] = df['attributes'].str.extract(r'gene_id=([^;]+)', expand=False)

        return df
    except Exception as e:
        print(f"Error loading or parsing GFF/GTF file {filepath}: {e}")
        return None

def image_to_base64(path):
    """Converts an image file to a base64 string for HTML embedding."""
    if not path or not Path(path).exists():
        return "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs=" # 1x1 transparent pixel
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
        return "N/A (No results file)"
    
    try:
        df = pd.read_csv(comp_file, sep='\t')
        top_hit = df[df['Query_ID'] == motif_id].sort_values('PCC', ascending=False).iloc[0]
        return f"{top_hit['DB_Name']} (PCC: {top_hit['PCC']:.2f})"
    except (IndexError, FileNotFoundError, KeyError):
        return "No match found"

def generate_html_report(data, output_path):
    """Generates the final HTML report file with visualization and table."""
    
    json_data = json.dumps(data, indent=4)
    
    db_headers = ""
    for db_name in data['db_names']:
        db_headers += f"<th>Top {db_name} Match</th>\n"

    table_rows_html = ""
    for motif in data['motifs_table']:
        db_matches_html = ""
        for db_name in data['db_names']:
            db_matches_html += f"<td>{motif['matches'].get(db_name, 'N/A')}</td>\n"
            
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
                </tr>
"""

    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>MOCA Report: {data['gene_id']}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 2em; }}
        h1, h2, h3 {{ color: #333; }}
        #chart-container {{ width: 100%; overflow-x: auto; border: 1px solid #ccc; background-color: #f9f9f9; margin-bottom: 2em;}}
        svg {{ min-width: 100%; }}
        .gene-backbone {{ stroke: #555; stroke-width: 2; }}
        .feature {{ cursor: pointer; }}
        .exon {{ fill: #B0B0B0; stroke: #888888; }} /* Light Grey for Exons */
        .CDS {{ fill: #666666; stroke: #444444; }} /* Dark Grey for CDS */
        .five_prime_UTR, .three_prime_UTR {{ fill: #E0E0E0; stroke: #C0C0C0; }} /* Very Light Grey for UTRs */
        .motif {{ cursor: pointer; }}
        .axis-line {{ stroke: #aaa; stroke-width: 1; }}
        .axis-tick {{ stroke: #aaa; stroke-width: 1; }}
        .axis-label {{ font-size: 10px; fill: #555; }}
        .tooltip {{ position: absolute; visibility: hidden; background: rgba(0, 0, 0, 0.8); color: #fff; padding: 8px; border-radius: 4px; font-size: 12px; pointer-events: none; white-space: pre; transition: opacity 0.2s; z-index: 10; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; font-size: 12px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; vertical-align: middle; }}
        th {{ background-color: #f2f2f2; }}
        img {{ max-width: 100%; height: auto; }}
        .logo-cell {{ width: 150px; }}
        .plot-cell {{ width: 250px; }}
    </style>
</head>
<body>
    <h1>MOCA Analysis Report</h1>
    <h2>Gene: <span id="gene-title"></span></h2>
    <h4>Chromosome: <span id="chr-title"></span> | Strand: <span id="gene-strand"></span> | Annotated Motifs: <span id="motif-count"></span></h4>
    
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
                <th>Genomic Position</th>
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
            document.getElementById('gene-title').textContent = data.gene_id;
            document.getElementById('chr-title').textContent = data.chr;
            document.getElementById('gene-strand').textContent = data.gene.strand;
            document.getElementById('motif-count').textContent = data.total_motifs;

            const margin = {{ top: 20, right: 50, bottom: 40, left: 50 }};
            const width = container.clientWidth - margin.left - margin.right;
            const height = 120 + (data.max_track * 30);

            const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
            svg.setAttribute("width", width + margin.left + margin.right);
            svg.setAttribute("height", height + margin.top + margin.bottom);
            
            const g = document.createElementNS("http://www.w3.org/2000/svg", "g");
            g.setAttribute("transform", `translate(${{margin.left}}, ${{margin.top}})`);
            svg.appendChild(g);
            container.appendChild(svg);

            const scale = (pos) => ((pos - data.view_start) / (data.view_end - data.view_start)) * width;

            // Draw axis
            const axisY = height - 20;
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
                label.setAttribute("y", axisY + 20);
                label.setAttribute("text-anchor", "middle");
                label.classList.add("axis-label");
                label.textContent = Math.round(pos); // Changed label to show bp instead of kb
                g.appendChild(label);
            }}

            const tooltip = document.getElementById('tooltip');

            // Draw gene model
            const geneY = 60;
            const backbone = document.createElementNS("http://www.w3.org/2000/svg", "line");
            backbone.setAttribute("x1", scale(data.gene.start));
            backbone.setAttribute("y1", geneY);
            backbone.setAttribute("x2", scale(data.gene.end));
            backbone.setAttribute("y2", geneY);
            backbone.classList.add("gene-backbone");
            g.appendChild(backbone);

            const numArrows = Math.floor((scale(data.gene.end) - scale(data.gene.start)) / 50);
            for (let i = 1; i < numArrows; i++) {{
                const arrowPos = scale(data.gene.start) + i * 50;
                const arrow = document.createElementNS("http://www.w3.org/2000/svg", "path");
                const arrowDir = data.gene.strand === '+' ? `M${{arrowPos-3}},${{geneY-3}} L${{arrowPos}},${{geneY}} L${{arrowPos-3}},${{geneY+3}}` : `M${{arrowPos+3}},${{geneY-3}} L${{arrowPos}},${{geneY}} L${{arrowPos+3}},${{geneY+3}}`;
                arrow.setAttribute("d", arrowDir);
                arrow.setAttribute("stroke", "#555");
                arrow.setAttribute("stroke-width", "1.5");
                arrow.setAttribute("fill", "none");
                g.appendChild(arrow);
            }}

            data.features.forEach(feature => {{
                const featureRect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
                featureRect.setAttribute("x", scale(feature.start));
                const featureWidth = scale(feature.end) - scale(feature.start);
                featureRect.setAttribute("width", Math.max(1, featureWidth));
                featureRect.classList.add("feature");
                featureRect.classList.add(feature.type);

                // Adjust height and Y position based on feature type
                if (feature.type === 'CDS') {{
                    featureRect.setAttribute("y", geneY - 7.5);
                    featureRect.setAttribute("height", 15);
                }} else {{ // UTRs and Exons
                    featureRect.setAttribute("y", geneY - 5);
                    featureRect.setAttribute("height", 10);
                }}
                g.appendChild(featureRect);

                featureRect.addEventListener('mouseover', (event) => {{
                    tooltip.style.visibility = 'visible';
                    tooltip.innerHTML = `<b>${{feature.type.replace('_', ' ')}}</b>\\n<b>Position:</b> ${{feature.start}} - ${{feature.end}}`;
                }});
                featureRect.addEventListener('mousemove', (event) => {{
                    tooltip.style.top = (event.pageY - 10) + 'px';
                    tooltip.style.left = (event.pageX + 10) + 'px';
                }});
                featureRect.addEventListener('mouseout', () => {{
                    tooltip.style.visibility = 'hidden';
                }});
            }});

            // Draw motifs
            data.motifs_viz.forEach(motif => {{
                const motifShape = document.createElementNS("http://www.w3.org/2000/svg", "path");
                const motifY = 85 + (motif.track * 20);
                const x1 = scale(motif.start);
                const x2 = scale(motif.end);
                
                let d;
                if (motif.strand === '+') {{
                    d = `M${{x1}},${{motifY}} L${{x2}},${{motifY + 5}} L${{x1}},${{motifY + 10}} Z`;
                }} else {{
                    d = `M${{x2}},${{motifY}} L${{x1}},${{motifY + 5}} L${{x2}},${{motifY + 10}} Z`;
                }}

                motifShape.setAttribute("d", d);
                motifShape.setAttribute("fill", motif.color);
                motifShape.classList.add("motif");
                g.appendChild(motifShape);

                motifShape.addEventListener('mouseover', (event) => {{
                    tooltip.style.visibility = 'visible';
                    tooltip.innerHTML = `<b>Motif:</b> ${{motif.id}}\\n<b>Score:</b> ${{motif.score.toFixed(2)}}\\n<b>Position:</b> ${{motif.start}} - ${{motif.end}}\\n<b>Strand:</b> ${{motif.strand}}`;
                }});
                motifShape.addEventListener('mousemove', (event) => {{
                    tooltip.style.top = (event.pageY - 10) + 'px';
                    tooltip.style.left = (event.pageX + 10) + 'px';
                }});
                motifShape.addEventListener('mouseout', () => {{
                    tooltip.style.visibility = 'hidden';
                }});
            }});
        }});
    </script>
</body>
</html>
    """
    with open(output_path, 'w') as f:
        f.write(html_template)

def run(config, common_settings):
    """
    Main function for the reporting step. Consolidates all analysis results
    for a target gene into a single HTML report.
    """
    # --- 1. Configuration & Path Setup ---
    output_dir = config['output_dir']
    
    annotation_dir = common_settings.get('annotation', {}).get('output_dir')
    viz_dir = common_settings.get('visualization', {}).get('output_dir')
    ranging_dir = common_settings.get('ranging', {}).get('output_dir')
    importance_dir = common_settings.get('importance', {}).get('output_dir')
    comparison_dir = common_settings.get('comparison', {}).get('output_dir')

    ref_gff_file = config.get('reference_gff')
    if not ref_gff_file:
        ref_gff_file = common_settings.get('annotation', {}).get('reference_gff')

    filter_method = common_settings.get('annotation', {}).get('filter_method', 'minmax')
    annotated_motifs_file = os.path.join(annotation_dir, f"annotated_motifs_{filter_method}.csv")
    importance_file = os.path.join(importance_dir, f"{common_settings.get('date', 'nodate')}_{common_settings.get('species_tag', 'unk')}{common_settings.get('model_tag', 'm0')}_contrib_scores.csv")

    if not ref_gff_file or not os.path.exists(ref_gff_file):
        print(f"Error: Reference GFF file not found at '{ref_gff_file}'.")
        return

    if not os.path.exists(annotated_motifs_file):
        print(f"Error: Annotated motifs file not found at '{annotated_motifs_file}'. Skipping report generation.")
        return

    # --- 2. Load Data and Select Target Gene(s) ---
    print("Loading annotated motifs for report...")
    annotated_df = pd.read_csv(annotated_motifs_file)
    if annotated_df.empty:
        print("Annotated motifs file is empty. Nothing to report.")
        return
        
    importance_df = pd.read_csv(importance_file) if os.path.exists(importance_file) else None
    
    target_gene_list_path = config.get('target_gene_list')
    genes_to_process = []
    if target_gene_list_path and os.path.exists(target_gene_list_path):
        print(f"Loading target genes from {target_gene_list_path}")
        with open(target_gene_list_path, 'r') as f:
            genes_to_process = [line.strip() for line in f if line.strip()]
    else:
        genes_to_process.append(annotated_df['gene_id'].iloc[0])

    print(f"Will generate reports for: {genes_to_process}")

    annotated_df['gene_id_norm'] = _normalize_gene_id(annotated_df['gene_id'])

    # --- Main loop to generate a report for each gene ---
    for target_gene_id_raw in genes_to_process:
        target_gene_id = re.sub(r'\.\d+$', '', target_gene_id_raw)
        target_gene_id = re.sub(r'^gene:', '', target_gene_id)
        print(f"\n--- Generating report for gene: {target_gene_id} ---")
        
        target_motifs_df = annotated_df[annotated_df['gene_id_norm'] == target_gene_id].copy()
        
        if target_motifs_df.empty:
            print(f"Warning: No annotated motifs found for gene '{target_gene_id}'. Skipping report.")
            continue
        
        # --- 3. Load Gene and Exon Information from GFF ---
        print(f"Loading features for {target_gene_id} from GFF...")
        gff_df = _load_gff(ref_gff_file)
        
        # Find the main gene entry using its normalized ID
        gff_df['gene_id_attr_norm'] = _normalize_gene_id(gff_df['gene_id_attr'])
        gene_info = gff_df[(gff_df['gene_id_attr_norm'] == target_gene_id) & (gff_df['type'] == 'gene')]
        
        if gene_info.empty:
            print(f"Error: Could not find gene '{target_gene_id}' in GFF file.")
            continue
            
        gene_gff_id = gene_info['id'].iloc[0]
        gene_start, gene_end = gene_info['start'].iloc[0], gene_info['end'].iloc[0]
        gene_chr = gene_info['seqid'].iloc[0]
        gene_strand = gene_info['strand'].iloc[0]
        
        # Find all transcripts whose Parent is the gene's GFF ID
        all_transcripts = gff_df[(gff_df['parent_id'] == gene_gff_id) & (gff_df['type'] == 'mRNA')]
        
        if all_transcripts.empty:
            print(f"Warning: No mRNA transcripts found for gene '{target_gene_id}'. No features will be displayed.")
            feature_info = pd.DataFrame()
        else:
            # Select the LONGEST transcript variant
            all_transcripts['length'] = all_transcripts['end'] - all_transcripts['start']
            longest_transcript = all_transcripts.loc[all_transcripts['length'].idxmax()]
            selected_transcript_id = longest_transcript['id']
            
            print(f"Found {len(all_transcripts)} transcript(s). Using the longest, '{selected_transcript_id}', for visualization.")
            
            # Find all features whose Parent is the selected transcript's ID, excluding CDS
            feature_types = ['exon', 'five_prime_UTR', 'three_prime_UTR']
            feature_info = gff_df[(gff_df['parent_id'] == selected_transcript_id) & (gff_df['type'].isin(feature_types))]
            print(f"Found {len(feature_info)} features (exons, UTRs) for this transcript.")

        # --- 4. Prepare Data for Visualization and Table ---
        target_motifs_df.sort_values('gen_mstart', inplace=True)
        tracks = []
        view_start = gene_start - 1500
        view_end = gene_end + 1500
        padding = (view_end - view_start) * 0.005

        for index, motif in target_motifs_df.iterrows():
            placed = False
            for i, track_end_pos in enumerate(tracks):
                if motif['gen_mstart'] > (track_end_pos + padding):
                    tracks[i] = motif['gen_mend']
                    target_motifs_df.loc[index, 'track'] = i
                    placed = True
                    break
            if not placed:
                tracks.append(motif['gen_mend'])
                target_motifs_df.loc[index, 'track'] = len(tracks) - 1

        unique_motifs = target_motifs_df['motif'].unique()
        colors = [f"hsl({(i * 360 / len(unique_motifs)) % 360}, 70%, 50%)" for i in range(len(unique_motifs))]
        color_map = dict(zip(unique_motifs, colors))
        target_motifs_df['color'] = target_motifs_df['motif'].map(color_map)

        # --- 5. Assemble data into a JSON-compatible dictionary ---
        db_names = [d.name for d in Path(comparison_dir).iterdir() if d.is_dir()] if os.path.exists(comparison_dir) else []
        if db_names:
            print(f"Found comparison results for databases: {db_names}")

        report_data = {
            'gene_id': target_gene_id,
            'chr': gene_chr,
            'view_start': int(view_start),
            'view_end': int(view_end),
            'gene': {'start': int(gene_start), 'end': int(gene_end), 'strand': gene_strand},
            'features': [{'start': int(r['start']), 'end': int(r['end']), 'type': r['type']} for _, r in feature_info.iterrows()],
            'motifs_viz': [{'id': r['motif'], 'start': int(r['gen_mstart']), 'end': int(r['gen_mend']), 
                            'score': r['score'], 'track': int(r['track']), 'color': r['color'], 'strand': r['strand']} 
                           for _, r in target_motifs_df.iterrows()],
            'max_track': int(target_motifs_df['track'].max()) if not target_motifs_df.empty and 'track' in target_motifs_df.columns else 0,
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
            
            logo_fwd_path = os.path.join(viz_dir, f"{base_motif_id}.png") if viz_dir else None
            logo_rev_path = os.path.join(viz_dir, f"{base_motif_id.replace('F_', 'R_')}.png") if viz_dir else None
            range_plot_path = os.path.join(ranging_dir, 'distribution_plots', f"epm_{common_settings.get('species_tag', 'unk')}_{common_settings.get('model_tag', 'm0')}_{motif_short_name}_density.png") if ranging_dir and motif_short_name else None

            importance_score = 0.0
            if importance_df is not None:
                score_row = importance_df[importance_df['motif'] == motif_id]
                if not score_row.empty:
                    importance_score = score_row['contrib_score_sum'].iloc[0]
            
            matches = {db_name: get_top_match(motif_id, comparison_dir, db_name) for db_name in db_names}
            
            if is_reverse:
                fwd_logo_for_display = image_to_base64(logo_rev_path)
                rev_logo_for_display = image_to_base64(logo_fwd_path)
            else:
                fwd_logo_for_display = image_to_base64(logo_fwd_path)
                rev_logo_for_display = image_to_base64(logo_rev_path)

            report_data['motifs_table'].append({
                'id': motif_id,
                'logo_fwd': fwd_logo_for_display,
                'logo_rev': rev_logo_for_display,
                'range_plot': image_to_base64(range_plot_path),
                'importance': importance_score,
                'position': f"{row['gen_mstart']}-{row['gen_mend']}",
                'distance': f"{row['dist_transc_border']}bp ({row['region']})",
                'strand': row['strand'],
                'matches': matches
            })

        # --- 6. Generate and Save HTML file ---
        output_path = os.path.join(output_dir, f"{target_gene_id}_moca_report.html")
        generate_html_report(report_data, output_path)

        print(f"Interactive HTML report saved to: {output_path}")
        
    print("\nReporting step successfully completed for all specified genes.")
