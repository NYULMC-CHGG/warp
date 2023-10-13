package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config.RNAWithUMIsConfig
import org.broadinstitute.dsp.pipelines.inputs.{
  RNAWithUMIsInputs,
  RNAWithUMIsValidationInputs
}

class RNAWithUMIsTester(testerConfig: RNAWithUMIsConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "RNAWithUMIsPipeline"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "rna_seq"

  override protected val validationWorkflowName: String =
    "VerifyRNAWithUMIs"

  // Validation uses the same options as the arrays workflow
  override protected lazy val validationWdlOptions: String = {
    readTestOptions(releaseDir, env)
  }

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/rna_seq/rna_with_umis/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/rna_seq/rna_with_umis/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val rnaWithUmisInputs = new RNAWithUMIsInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName =
      rnaWithUmisInputs.getBasename(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath

    val validationInputs = RNAWithUMIsValidationInputs(
      test_metrics = Array(
        resultsCloudPath.resolve(
          s"$outputBaseName.transcriptome.duplicate.metrics"),
        resultsCloudPath.resolve(s"$outputBaseName.alignment_summary_metrics"),
        resultsCloudPath.resolve(s"$outputBaseName.duplicate.metrics"),
        resultsCloudPath.resolve(s"$outputBaseName.insert_size_metrics"),
        resultsCloudPath.resolve(s"$outputBaseName.rna_metrics"),
        resultsCloudPath.resolve(
          s"$outputBaseName.quality_distribution_metrics"),
        resultsCloudPath.resolve(s"$outputBaseName.quality_by_cycle_metrics"),
        resultsCloudPath.resolve(
          s"$outputBaseName.base_distribution_by_cycle_metrics")
      ),
      truth_metrics = Array(
        truthCloudPath.resolve(
          s"$outputBaseName.transcriptome.duplicate.metrics"),
        truthCloudPath.resolve(s"$outputBaseName.alignment_summary_metrics"),
        truthCloudPath.resolve(s"$outputBaseName.duplicate.metrics"),
        truthCloudPath.resolve(s"$outputBaseName.insert_size_metrics"),
        truthCloudPath.resolve(s"$outputBaseName.rna_metrics"),
        truthCloudPath.resolve(s"$outputBaseName.quality_distribution_metrics"),
        truthCloudPath.resolve(s"$outputBaseName.quality_by_cycle_metrics"),
        truthCloudPath.resolve(
          s"$outputBaseName.base_distribution_by_cycle_metrics")
      ),
      test_text_metrics = Array(
        resultsCloudPath.resolve(s"$outputBaseName.metrics.tsv")
      ),
      truth_text_metrics = Array(
        truthCloudPath.resolve(s"$outputBaseName.metrics.tsv")
      ),
      test_output_bam = resultsCloudPath.resolve(
        s"$outputBaseName.duplicate_marked.coordinate_sorted.bam"),
      truth_output_bam = truthCloudPath.resolve(
        s"$outputBaseName.duplicate_marked.coordinate_sorted.bam"),
      test_transcriptome_bam = resultsCloudPath.resolve(
        s"$outputBaseName.transcriptome_RSEM_post_processed.bam"),
      truth_transcriptome_bam = truthCloudPath.resolve(
        s"$outputBaseName.transcriptome_RSEM_post_processed.bam"),
      test_gene_tpm =
        resultsCloudPath.resolve(s"$outputBaseName.gene_tpm.gct.gz"),
      truth_gene_tpm =
        truthCloudPath.resolve(s"$outputBaseName.gene_tpm.gct.gz"),
      test_gene_counts =
        resultsCloudPath.resolve(s"$outputBaseName.gene_reads.gct.gz"),
      truth_gene_counts =
        truthCloudPath.resolve(s"$outputBaseName.gene_reads.gct.gz"),
      test_exon_counts =
        resultsCloudPath.resolve(s"$outputBaseName.exon_reads.gct.gz"),
      truth_exon_counts =
        truthCloudPath.resolve(s"$outputBaseName.exon_reads.gct.gz")
    )
    RNAWithUMIsValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
