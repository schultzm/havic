def snakefile(**kwargs):
    """
    This will generate the Snakefile for running with snakemake
    :param kwargs:
    :return: Snakefile
    """
    output:


    # Pipeline starts here with Ruffus
    @mkdir(self.outdir)
    def create_outdir():
        pass
        # print(f"Creating output directory {self.outdir}")


    @follows(create_outdir)
    @files(self.query_files,
           self.outfiles['tmp_fasta'], self.havnet_ampliconseq)
    def compile_input_fasta(infile, outfile, refamplicon):
        self._compile_input_fasta()

    
    @follows(compile_input_fasta)
    @files(self.outfiles['tmp_fasta'],
           self.outfiles['tmp_bam'])
    def minimap2_input_fasta_to_ref(infile, outfile):
        self._minimap2_input_fasta_to_ref()


    @follows(minimap2_input_fasta_to_ref)
    @files(self.outfiles['tmp_bam'],
           self.outfiles['fasta_from_bam'])
    def bam2fasta(infile, outfile):
        self._bam2fasta()


    @follows(bam2fasta)
    @files(self.outfiles['fasta_from_bam'],
           self.outfiles['fasta_from_bam_trimmed'])
    def get_cleaned_fasta(infile, outfile):
        aln = self._get_clean_fasta_alignment()
        if aln and len(aln) < 3:
            exit_statement = f'{aln}\n' + \
                             f'Need at least three sequences in ' + \
                             f'alignment to continue (n={len(aln)})'
            import sys
            sys.exit(exit_statement)


    @follows(get_cleaned_fasta)
    @files(self.outfiles['fasta_from_bam_trimmed'],
           self.outfiles["mp_treefile"])
    def run_iqtree(infile, outfile):
        self._run_iqtree()


    @follows(run_iqtree)
    @files(self.outfiles['treefile'],
           self.outfiles['mp_treefile'])
    def midpoint_root_iqtree(infile, outfile):
        self._midpoint_root_iqtree()


    @follows(midpoint_root_iqtree)
    @files([self.outfiles['fasta_from_bam_trimmed'],
            self.outfiles['mp_treefile']],
           self.outfiles["clusterpicked_tree"])
    def clusterpick_from_mpr_iqtree_and_cleaned_fasta(infile, outfile):
        self._clusterpick()


    @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
    @files(self.outfiles["cluster_assignments"],
           self.outfiles["clusters_assigned"])
    def summarise_cluster_assignments(infile, outfile):
        self._summarise_cluster_assignments()


    @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
    @files(self.outfiles['clusterpicked_tree'],
           self.outfiles['clusterpicked_tree_bkp'])
    def backup_clusterpicked_figtree(infile, outfile):
        self._bkp_clusterpickedtree()


    @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
    @files([self.outfiles['fasta_from_bam_trimmed'],
            self.outfiles['mp_treefile']],
           self.outfiles['treeplotr'])
    def plot_results_ggtree(infiles, outfiles):
        self._plot_results()


    # Run the pipeline
    import tempfile


    def mv_tmp_sqlite(temp_sqlite_db, perm_sqlite_db):
        try:
            os.popen(f'mv {temp_sqlite_db} {perm_sqlite_db}')
            print(f"Moved {temp_sqlite_db} to {perm_sqlite_db}.")
        except IOError:
            print(f'Unable to move {temp_sqlite_db} to {self.outdir}',
                  file=sys.stderr)


    with tempfile.TemporaryDirectory() as tmpfile:
        # final_outfolder = self.outdir
        # self.outdir = tmpfile
        db_name = '.ruffus_history.sqlite'
        temp_sqlite_db = os.path.join(tmpfile, db_name)
        perm_sqlite_db = os.path.join(os.path.abspath(self.outdir),
                                      db_name)
        if os.path.exists(perm_sqlite_db):
            os.popen(f'cp {perm_sqlite_db} {temp_sqlite_db}')
            print(f'Copied {perm_sqlite_db} to {temp_sqlite_db}.')
        else:
            print(f'Making new SQLite db at {temp_sqlite_db}')
        if self.redo:
            pipeline_run(forcedtorun_tasks=compile_input_fasta,
                         history_file=temp_sqlite_db)
            mv_tmp_sqlite(temp_sqlite_db, perm_sqlite_db)
        else:
            pipeline_run(history_file=temp_sqlite_db)
            mv_tmp_sqlite(temp_sqlite_db, perm_sqlite_db)

        # Print out the pipeline graph
        pipeline_printout_graph(
            make_path(self.outdir, self.prefix, "_pipeline_graph.svg"),
            "svg")