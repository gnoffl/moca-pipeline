import os
import pandas as pd
import json
import re
import numpy as np
import base64
from pathlib import Path

def _load_gff(filepath):
    """Loads a GFF/GTF file, parsing it into a pandas DataFrame."""
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        
        df['gene_id'] = df['attributes'].str.extract(r'gene_id[= ]"([^"]+)"', expand=False)
        if df['gene_id'].isnull().all():
            df['gene_id'] = df['attributes'].str.extract(r'ID=([^;]+)', expand=False)
        
        df.dropna(subset=['gene_id'], inplace=True)
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
    
    # Dynamically create table headers from the discovered databases
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
        .exon {{ fill: #3498db; stroke: #2980b9; cursor: pointer; }}
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
                label.textContent = Math.round(pos / 1000) + 'kb';
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

            data.exons.forEach(exon => {{
                const exonRect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
                exonRect.setAttribute("x", scale(exon.start));
                exonRect.setAttribute("y", geneY - 7.5);
                exonRect.setAttribute("width", scale(exon.end) - scale(exon.start));
                exonRect.setAttribute("height", 15);
                exonRect.classList.add("exon");
                g.appendChild(exonRect);

                exonRect.addEventListener('mouseover', (event) => {{
                    tooltip.style.visibility = 'visible';
                    tooltip.innerHTML = `<b>Exon</b>\\n<b>Position:</b> ${{exon.start}} - ${{exon.end}}`;
                }});
                exonRect.addEventListener('mousemove', (event) => {{
                    tooltip.style.top = (event.pageY - 10) + 'px';
                    tooltip.style.left = (event.pageX + 10) + 'px';
                }});
                exonRect.addEventListener('mouseout', () => {{
                    tooltip.style.visibility = 'hidden';
                }});
            }});

            // Draw motifs
            data.motifs_viz.forEach(motif => {{
                const motifShape = document.createElementNS("http://www.w3.org/2000/svg", "path");
                const motifY = 85 + (motif.track * 20);
                const x1 = scale(motif.start);
                const x2 = scale(motif.end);
                const w = Math.max(1, x2 - x1);
                
                let d;
                if (w < 8) {{
                    d = `M${{x1}},${{motifY}} h${{w}} v10 h-${{w}} z`;
                }} else {{
                    if (motif.strand === '+') {{
                        d = `M${{x1}},${{motifY}} h${{w-5}} l5,5 l-5,5 h-${{w-5}} z`;
                    }} else {{
                        d = `M${{x2}},${{motifY}} h-${{w-5}} l-5,5 l5,5 h${{w-5}} z`;
                    }}
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
    ref_gff_file = common_settings.get('annotation', {}).get('reference_gff')
    
    filter_method = common_settings.get('annotation', {}).get('filter_method', 'q1q9')
    annotated_motifs_file = os.path.join(annotation_dir, f"annotated_motifs_{filter_method}.csv")
    importance_file = os.path.join(importance_dir, f"{common_settings.get('date', 'nodate')}_{common_settings.get('species_tag', 'unk')}{common_settings.get('model_tag', 'm0')}_contrib_scores.csv")

    if not os.path.exists(annotated_motifs_file):
        print("Error: Annotated motifs file not found. Skipping report generation.")
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
        # Default to the first gene if no list is provided
        genes_to_process.append(annotated_df['gene_id'].iloc[0])

    print(f"Will generate reports for: {genes_to_process}")

    # --- Main loop to generate a report for each gene ---
    for target_gene_id in genes_to_process:
        print(f"\n--- Generating report for gene: {target_gene_id} ---")
        
        target_motifs_df = annotated_df[annotated_df['gene_id'] == target_gene_id].copy()
        if target_motifs_df.empty:
            print(f"Warning: No annotated motifs found for gene '{target_gene_id}'. Skipping report.")
            continue
        
        # --- 3. Load Gene and Exon Information from GFF ---
        print(f"Loading features for {target_gene_id} from GFF...")
        gff_df = _load_gff(ref_gff_file)
        
        gene_info = gff_df[(gff_df['gene_id'] == target_gene_id) & (gff_df['type'] == 'gene')]
        if gene_info.empty:
            print(f"Error: Could not find gene '{target_gene_id}' in GFF file.")
            continue
            
        gene_start, gene_end = gene_info['start'].iloc[0], gene_info['end'].iloc[0]
        gene_chr = gene_info['seqid'].iloc[0]
        gene_strand = gene_info['strand'].iloc[0]
        
        exon_info = gff_df[(gff_df['gene_id'] == target_gene_id) & (gff_df['type'] == 'exon')]

        # --- 4. Prepare Data for Visualization and Table ---
        target_motifs_df.sort_values('gen_mstart', inplace=True)
        tracks = []
        view_start = gene_start - 1000
        view_end = gene_end + 1000
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
        db_names = [d.name for d in Path(comparison_dir).iterdir() if d.is_dir()]
        print(f"Found comparison results for databases: {db_names}")

        report_data = {
            'gene_id': target_gene_id,
            'chr': gene_chr,
            'view_start': int(view_start),
            'view_end': int(view_end),
            'gene': {'start': int(gene_start), 'end': int(gene_end), 'strand': gene_strand},
            'exons': [{'start': int(r['start']), 'end': int(r['end'])} for _, r in exon_info.iterrows()],
            'motifs_viz': [{'id': r['motif'], 'start': int(r['gen_mstart']), 'end': int(r['gen_mend']), 
                            'score': r['score'], 'track': int(r['track']), 'color': r['color'], 'strand': r['strand']} 
                           for _, r in target_motifs_df.iterrows()],
            'max_track': int(target_motifs_df['track'].max()) if not target_motifs_df.empty else 0,
            'total_motifs': len(target_motifs_df),
            'db_names': db_names,
            'motifs_table': []
        }

        for _, row in target_motifs_df.iterrows():
            motif_id = row['motif']
            is_reverse = 'R_' in motif_id
            base_motif_id = motif_id.replace('R_', 'F_')
            motif_short_name = re.search(r'(p\d+m\d+)', base_motif_id).group(1) if re.search(r'(p\d+m\d+)', base_motif_id) else ""
            
            logo_fwd_path = os.path.join(viz_dir, f"{base_motif_id}.png")
            logo_rev_path = os.path.join(viz_dir, f"{base_motif_id.replace('F_', 'R_')}.png")
            range_plot_path = os.path.join(ranging_dir, 'distribution_plots', f"epm_{common_settings['species_tag']}_{common_settings['model_tag']}_{motif_short_name}_density.png")
            
            importance_score = 0.0
            if importance_df is not None:
                score_row = importance_df[importance_df['motif'] == motif_id]
                if not score_row.empty:
                    importance_score = score_row['contrib_score_sum'].iloc[0]
            
            matches = {db_name: get_top_match(motif_id, comparison_dir, db_name) for db_name in db_names}
            
            # Swap the logos if the motif is a reverse match
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
