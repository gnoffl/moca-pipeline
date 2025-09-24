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

def _get_deepcre_pos(genomic_pos, tss, tts, gene_strand):
    """Calculates the deepCRE coordinate (1-3000) for a given genomic position."""
    try:
        genomic_pos, tss, tts = int(genomic_pos), int(tss), int(tts)
    except (ValueError, TypeError):
        return None

    if gene_strand == '+':
        if tss - 1000 <= genomic_pos < tss:
            return (genomic_pos - (tss - 1000)) + 1
        elif tss <= genomic_pos < tss + 500:
            return (genomic_pos - tss) + 1001
        elif tts - 500 <= genomic_pos < tts:
            return (genomic_pos - (tts - 500)) + 1501
        elif tts <= genomic_pos < tts + 1000:
            return (genomic_pos - tts) + 2001
    else:  # gene_strand == '-'
        if tss < genomic_pos <= tss + 1000:
            return (tss + 1000 - genomic_pos) + 1
        elif tss - 500 < genomic_pos <= tss:
            return (tss - genomic_pos) + 1001
        elif tts < genomic_pos <= tts + 500:
            return (tts + 500 - genomic_pos) + 1501
        elif tts - 1000 < genomic_pos <= tts:
            return (tts - genomic_pos) + 2001
    return None

def _load_gff(filepath):
    """Loads a GFF/GTF file, parsing it into a pandas DataFrame."""
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        
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
                    <td>{motif['deepcre_pos']}</td>
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
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 2em; background-color: #fdfdfd; }}
        h1, h2, h3 {{ color: #333; }}
        #chart-container {{ width: 100%; overflow-x: auto; border: 1px solid #ccc; background-color: #f9f9f9; margin-bottom: 2em; }}
        svg {{ min-width: 100%; }}
        .gene-backbone {{ stroke: #555; stroke-width: 2; }}
        .feature {{ cursor: pointer; }}
        .exon {{ fill: #666666; stroke: #444444; }}
        .CDS {{ fill: #666666; stroke: #444444; }}
        .five_prime_UTR, .three_prime_UTR {{ fill: #a0a0a0; stroke: #888888; }}
        .motif-group {{ cursor: pointer; }}
        .axis-line, .axis-tick {{ stroke: #aaa; stroke-width: 1; }}
        .axis-label {{ font-size: 12px; fill: #555; }}
        .deepcre-axis-segment {{ stroke: #333; stroke-width: 1.5; }}
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
                <th>deepCRE Position</th>
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

            const margin = {{ top: 60, right: 50, bottom: 60, left: 50 }};
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
            const backbone = document.createElementNS("http://www.w3.org/2000/svg", "line");
            backbone.setAttribute("x1", scale(data.gene.start));
            backbone.setAttribute("y1", geneY);
            backbone.setAttribute("x2", scale(data.gene.end));
            backbone.setAttribute("y2", geneY);
            backbone.classList.add("gene-backbone");
            g.appendChild(backbone);
            
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
                featureRect.addEventListener('mouseover', (event) => {{
                    tooltip.style.visibility = 'visible';
                    tooltip.innerHTML = `<b>${{feature.type.replace('_', ' ')}}</b>\\n<b>Position:</b> ${{feature.start.toLocaleString()}} - ${{feature.end.toLocaleString()}}`;
                }});
                featureRect.addEventListener('mousemove', (event) => {{
                    tooltip.style.top = (event.pageY - 10) + 'px';
                    tooltip.style.left = (event.pageX + 10) + 'px';
                }});
                featureRect.addEventListener('mouseout', () => {{
                    tooltip.style.visibility = 'hidden';
                }});
            }});
            
            data.gene_labels.forEach(label => {{
                const labelText = document.createElementNS("http://www.w3.org/2000/svg", "text");
                labelText.setAttribute("x", scale(label.pos));
                labelText.setAttribute("y", geneY - 25);
                labelText.setAttribute("text-anchor", "middle");
                labelText.setAttribute("font-size", "12px");
                labelText.setAttribute("font-weight", "bold");
                labelText.textContent = label.text;
                g.appendChild(labelText);
            }});

            // Draw motifs
            data.motifs_viz.forEach(motif => {{
                const motifGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
                motifGroup.classList.add("motif-group");
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
                    tooltip.innerHTML = `<b>Motif:</b> ${{motif.id}}\\n<b>Score:</b> ${{motif.score.toFixed(2)}}\\n<b>Position:</b> ${{motif.start.toLocaleString()}} - ${{motif.end.toLocaleString()}}\\n<b>Strand:</b> ${{motif.strand}}`;
                }});
                motifGroup.addEventListener('mousemove', (event) => {{
                    tooltip.style.top = (event.pageY - 10) + 'px';
                    tooltip.style.left = (event.pageX + 10) + 'px';
                }});
                motifGroup.addEventListener('mouseout', () => {{
                    tooltip.style.visibility = 'hidden';
                }});
            }});

            // --- Draw Top Genomic Axis ---
            const topAxisY = -30;
            const topAxisLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
            topAxisLine.setAttribute("x1", 0);
            topAxisLine.setAttribute("y1", topAxisY);
            topAxisLine.setAttribute("x2", width);
            topAxisLine.setAttribute("y2", topAxisY);
            topAxisLine.classList.add("axis-line");
            g.appendChild(topAxisLine);

            const viewRange = data.view_end - data.view_start;
            const topTickInterval = Math.pow(10, Math.floor(Math.log10(viewRange)) - 1) * 5;
            let topCurrentPos = Math.ceil(data.view_start / topTickInterval) * topTickInterval;
            
            while (topCurrentPos <= data.view_end) {{
                const x = scale(topCurrentPos);
                const tick = document.createElementNS("http://www.w3.org/2000/svg", "line");
                tick.setAttribute("x1", x);
                tick.setAttribute("y1", topAxisY);
                tick.setAttribute("x2", x);
                tick.setAttribute("y2", topAxisY - 5);
                tick.classList.add("axis-tick");
                g.appendChild(tick);

                const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
                label.setAttribute("x", x);
                label.setAttribute("y", topAxisY - 10);
                label.setAttribute("text-anchor", "middle");
                label.classList.add("axis-label");
                label.textContent = topCurrentPos.toLocaleString();
                g.appendChild(label);
                
                topCurrentPos += topTickInterval;
            }}

            // --- Draw Bottom deepCRE Axis ---
            const bottomAxisY = height;
            const tss = data.tss_pos;
            const tts = data.tts_pos;
            const geneStrand = data.gene.strand;

            const keyPoints = [];
            if (geneStrand === '+') {{
                keyPoints.push({{ pos: tss - 1000, label: "1" }});
                keyPoints.push({{ pos: tss, label: "1001" }});
                keyPoints.push({{ pos: tss + 500, label: "1500" }});
                keyPoints.push({{ pos: tts - 500, label: "1501" }});
                keyPoints.push({{ pos: tts, label: "2001" }});
                keyPoints.push({{ pos: tts + 1000, label: "3000" }});
            }} else {{ // strand is '-'
                keyPoints.push({{ pos: tss + 1000, label: "1" }});
                keyPoints.push({{ pos: tss, label: "1001" }});
                keyPoints.push({{ pos: tss - 500, label: "1500" }});
                keyPoints.push({{ pos: tts + 500, label: "1501" }});
                keyPoints.push({{ pos: tts, label: "2001" }});
                keyPoints.push({{ pos: tts - 1000, label: "3000" }});
            }}

            const drawSegment = (startPos, endPos) => {{
                const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
                line.setAttribute("x1", scale(startPos));
                line.setAttribute("y1", bottomAxisY);
                line.setAttribute("x2", scale(endPos));
                line.setAttribute("y2", bottomAxisY);
                line.classList.add("deepcre-axis-segment");
                g.appendChild(line);
            }};

            drawSegment(keyPoints[0].pos, keyPoints[1].pos); // Upstream
            drawSegment(keyPoints[1].pos, keyPoints[2].pos); // Gene Start
            drawSegment(keyPoints[3].pos, keyPoints[4].pos); // Gene End
            drawSegment(keyPoints[4].pos, keyPoints[5].pos); // Downstream

            keyPoints.forEach(point => {{
                if (point.pos >= data.view_start && point.pos <= data.view_end) {{
                    const x = scale(point.pos);
                    const tick = document.createElementNS("http://www.w3.org/2000/svg", "line");
                    tick.setAttribute("x1", x);
                    tick.setAttribute("y1", bottomAxisY);
                    tick.setAttribute("x2", x);
                    tick.setAttribute("y2", bottomAxisY + 5);
                    tick.classList.add("axis-tick");
                    g.appendChild(tick);
                    
                    const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
                    label.setAttribute("x", x);
                    label.setAttribute("y", bottomAxisY + 20);
                    label.setAttribute("text-anchor", "middle");
                    label.classList.add("axis-label");
                    label.textContent = point.label;
                    g.appendChild(label);
                }}
            }});
        }});
    </script>
</body>
</html>
    """
    with open(output_path, 'w', encoding='utf-8') as f:
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
        
        gff_df['gene_id_attr_norm'] = _normalize_gene_id(gff_df['gene_id_attr'])
        gene_info = gff_df[(gff_df['gene_id_attr_norm'] == target_gene_id) & (gff_df['type'] == 'gene')]
        
        if gene_info.empty:
            print(f"Error: Could not find gene '{target_gene_id}' in GFF file.")
            continue
            
        gene_gff_id = gene_info['id'].iloc[0]
        gene_start, gene_end = gene_info['start'].iloc[0], gene_info['end'].iloc[0]
        gene_chr = gene_info['seqid'].iloc[0]
        gene_strand = gene_info['strand'].iloc[0]
        
        # Define TSS (Transcription Start Site) and TTS (Transcription Termination Site)
        if gene_strand == '+':
            tss_pos = gene_start
            tts_pos = gene_end
        else:
            tss_pos = gene_end
            tts_pos = gene_start
        
        all_transcripts = gff_df[(gff_df['parent_id'] == gene_gff_id) & (gff_df['type'] == 'mRNA')]
        
        if all_transcripts.empty:
            print(f"Warning: No mRNA transcripts found for gene '{target_gene_id}'. No features will be displayed.")
            feature_info = pd.DataFrame()
        else:
            all_transcripts['length'] = all_transcripts['end'] - all_transcripts['start']
            longest_transcript = all_transcripts.loc[all_transcripts['length'].idxmax()]
            selected_transcript_id = longest_transcript['id']
            
            print(f"Found {len(all_transcripts)} transcript(s). Using the longest, '{selected_transcript_id}', for visualization.")
            
            feature_types = ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
            feature_info = gff_df[(gff_df['parent_id'] == selected_transcript_id) & (gff_df['type'].isin(feature_types))]
            print(f"Found {len(feature_info)} features (Exons, CDS, UTRs) for this transcript.")

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

        print("Applying custom color scheme to motifs...")
        unique_motifs = target_motifs_df['motif'].unique()
        color_map = {}
        
        p0_motifs = sorted([m for m in unique_motifs if '_p0' in m])
        p1_motifs = sorted([m for m in unique_motifs if '_p1' in m])
        other_motifs = sorted([m for m in unique_motifs if '_p0' not in m and '_p1' not in m])
        
        if p0_motifs:
            p0_hues = np.linspace(0, 60, len(p0_motifs), endpoint=True)
            for i, motif in enumerate(p0_motifs):
                color_map[motif] = f"hsl({p0_hues[i]}, 90%, 50%)"
                
        if p1_motifs:
            p1_hues = np.linspace(195, 255, len(p1_motifs), endpoint=True)
            for i, motif in enumerate(p1_motifs):
                color_map[motif] = f"hsl({p1_hues[i]}, 85%, 55%)"

        if other_motifs:
            other_lightness = np.linspace(40, 70, len(other_motifs), endpoint=True)
            for i, motif in enumerate(other_motifs):
                color_map[motif] = f"hsl(0, 0%, {other_lightness[i]:.0f}%)"
        
        target_motifs_df['color'] = target_motifs_df['motif'].map(color_map)

        print(f"Filtering motifs for visualization to match gene strand ('{gene_strand}')...")
        viz_motifs_df = target_motifs_df[target_motifs_df['strand'] == gene_strand].copy()
        print(f"Visualizing {len(viz_motifs_df)} of {len(target_motifs_df)} total motifs.")

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
            'tss_pos': int(tss_pos),
            'tts_pos': int(tts_pos),
            'features': [{'start': int(r['start']), 'end': int(r['end']), 'type': r['type']} for _, r in feature_info.iterrows()],
            'gene_labels': [{'pos': (gene_start + gene_end) / 2, 'text': target_gene_id}],
            'motifs_viz': [{'id': r['motif'], 'start': int(r['gen_mstart']), 'end': int(r['gen_mend']), 
                            'score': r['score'], 'track': int(r['track']), 'color': r['color'], 'strand': r['strand']} 
                           for _, r in viz_motifs_df.iterrows()],
            'max_track': int(viz_motifs_df['track'].max()) if not viz_motifs_df.empty and 'track' in viz_motifs_df.columns else 0,
            'total_motifs': len(target_motifs_df),
            'db_names': db_names,
            'motifs_table': []
        }

        # Create the detailed table data
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
            
            fwd_logo_for_display = image_to_base64(logo_rev_path if is_reverse else logo_fwd_path)
            rev_logo_for_display = image_to_base64(logo_fwd_path if is_reverse else logo_rev_path)

            motif_midpoint = (row['gen_mstart'] + row['gen_mend']) / 2
            deepcre_pos_val = _get_deepcre_pos(motif_midpoint, tss_pos, tts_pos, gene_strand)
            deepcre_pos_str = f"{deepcre_pos_val:.0f}" if deepcre_pos_val is not None else "N/A"

            report_data['motifs_table'].append({
                'id': motif_id,
                'logo_fwd': fwd_logo_for_display,
                'logo_rev': rev_logo_for_display,
                'range_plot': image_to_base64(range_plot_path),
                'importance': importance_score,
                'position': f"{row['gen_mstart']}-{row['gen_mend']}",
                'deepcre_pos': deepcre_pos_str,
                'strand': row['strand'],
                'matches': matches
            })

        # --- 6. Generate and Save HTML file ---
        output_path = os.path.join(output_dir, f"{target_gene_id}_moca_report.html")
        generate_html_report(report_data, output_path)

        print(f"Interactive HTML report saved to: {output_path}")
        
    print("\nReporting step successfully completed for all specified genes.")
