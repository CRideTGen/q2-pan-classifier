from os import path

from jinja2 import Environment, FileSystemLoader

from q2_pan_classifier.utils import merge_table, get_cpus


def classify_reads(
    ctx,
    trained_classifier,
    samp_reads=None,
    trunc_len_f=None,
    trunc_len_r=None,
    dada2_table=None,
    dada2_rep_seqs=None,
    dada2_stats=None,
):
    results = []

    # action importing
    dada2 = ctx.get_action("dada2", "denoise_paired")
    classify_sklearn = ctx.get_action("feature_classifier", "classify_sklearn")
    barplot = ctx.get_action("taxa", "barplot")
    transpose = ctx.get_action("feature_table", "transpose")
    tabulate = ctx.get_action("metadata", "tabulate")
    vis_test = ctx.get_action("pan_classifier", "visualization_final")

    # getting some output

    if dada2_table and dada2_rep_seqs and dada2_stats:
        dada2_table_out = dada2_table
        dada2_rep_seqs_out = dada2_rep_seqs
        dada2_stats_out = dada2_stats
    else:
        dada2_table_out, dada2_rep_seqs_out, dada2_stats_out = dada2(
            demultiplexed_seqs=samp_reads,
            trunc_len_f=trunc_len_f,
            trunc_len_r=trunc_len_r,
            n_threads=get_cpus(),
        )

    (classified,) = classify_sklearn(
        classifier=trained_classifier, reads=dada2_rep_seqs_out
    )
    barplot_taxonomy = barplot(table=dada2_table_out, taxonomy=classified)

    merged_table = merge_table(
        transpose, dada2_table_out, dada2_rep_seqs_out, classified
    )
    tabulated_table = tabulate(merged_table)

    results += [classified]
    results += barplot_taxonomy
    results += tabulated_table
    results += [dada2_table_out, dada2_rep_seqs_out, dada2_stats_out]

    return tuple(results)


def classify_reads_single(
    ctx,
    trained_classifier,
    samp_reads=None,
    trunc_len=None,
    dada2_table=None,
    dada2_rep_seqs=None,
    dada2_stats=None,
):
    results = []

    # action importing
    dada2 = ctx.get_action("dada2", "denoise_single")
    classify_sklearn = ctx.get_action("feature_classifier", "classify_sklearn")
    barplot = ctx.get_action("taxa", "barplot")
    transpose = ctx.get_action("feature_table", "transpose")
    tabulate = ctx.get_action("metadata", "tabulate")
    vis_test = ctx.get_action("pan_classifier", "visualization_final")

    # getting some output

    if dada2_table and dada2_rep_seqs and dada2_stats:
        dada2_table_out = dada2_table
        dada2_rep_seqs_out = dada2_rep_seqs
        dada2_stats_out = dada2_stats
    else:
        dada2_table_out, dada2_rep_seqs_out, dada2_stats_out = dada2(
            demultiplexed_seqs=samp_reads, trunc_len=trunc_len, n_threads=get_cpus()
        )

    (classified,) = classify_sklearn(
        classifier=trained_classifier, reads=dada2_rep_seqs_out
    )
    barplot_taxonomy = barplot(table=dada2_table_out, taxonomy=classified)

    merged_table = merge_table(
        transpose, dada2_table_out, dada2_rep_seqs_out, classified
    )
    tabulated_table = tabulate(merged_table)

    results += [classified]
    results += barplot_taxonomy
    results += tabulated_table
    results += [dada2_table_out, dada2_rep_seqs_out, dada2_stats_out]

    return tuple(results)


def visualization_final(output_dir: str) -> None:
    # temp_dir = tempfile.TemporaryDirectory()
    template_data = pkg_resources.resource_filename("q2_pan_classifier", "templates")
    jin_env = Environment(loader=FileSystemLoader(template_data), auto_reload=True)
    # shutil.copy2(path.join(template_data, 'base.html'), temp_dir.name)

    # jin_env = Environment(loader=FileSystemLoader("templates"))
    jin_out = jin_env.get_template("base.html").render(title="This is my title")

    with open(path.join(output_dir, "index.html"), "w") as f:
        f.write(jin_out)



