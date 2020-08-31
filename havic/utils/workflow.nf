#!/usr/bin/env nextflow


process createPipeline {
    
}
detection_pipeline = Pipeline(self.yaml)
# for key, value in detection_pipeline.__dict__.items():
#     print(f"{key}: {value}\n")
self.assertEqual(detection_pipeline.yaml_in['QUERY_FILES'][0],
                    'data/example1.fa')
detection_pipeline._run()

    def _run(self):
        """
        Run the pipeline using Ruffus.

        :return: None
        """

        os.system(f"nextflow run {rf(__parent_dir__, __nextflow_nf__)}")
        # Pipeline starts here without Ruffus
        # For development and testing.
        # os.mkdir(self.outdir)
        # self._compile_input_fasta()
        # self._minimap2_input_fasta_to_ref()
        # self._bam2fasta()
        # self._get_clean_fasta_alignment()
        # self._run_iqtree()
        # self._midpoint_root_iqtree()
        # self._clusterpick()
        # self._plot_results()

        # Pipeline starts here with Ruffus
        @mkdir(self.outdir)
        def create_outdir():
            pass

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
                exit_statement = f'{aln}\n' +\
                                 f'Need at least three sequences in ' +\
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
        with tempfile.TemporaryDirectory() as tmpfile:
            db_name = '.ruffus_history.sqlite'
            temp_sqlite_db = Path(tmpfile).joinpath(db_name)
            perm_sqlite_db = Path(self.outdir).joinpath(db_name)
            if self.yaml_in['FORCE_OVERWRITE_AND_RE_RUN']:
                pipeline_run(forcedtorun_tasks=create_outdir,
                             history_file=temp_sqlite_db)
                shutil.copyfile(temp_sqlite_db, perm_sqlite_db)
            else:
                shutil.copyfile(perm_sqlite_db, temp_sqlite_db)
                pipeline_run(history_file=temp_sqlite_db)
                shutil.copyfile(temp_sqlite_db, perm_sqlite_db)

            # Print out the pipeline graph
            pipeline_printout_graph(
                make_path(self.outdir, "pipeline_graph.svg"),
                "svg")