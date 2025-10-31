# JUBE with TSMP

JUBE (https://apps.fz-juelich.de/jsc/jube/jube2/docu/) is a JSC-developed tool to set up and analyze reproducible benchmarks. To set it up, JUBE must be installed on the system and configured with an *XML* or *YAML* file.

HPSCTerrSys has configured an XML file that enables JUBE-benchmarks with the available test cases of TSMP. Currently, the file is located at https://gitlab.jsc.fz-juelich.de/sdlts/benchmarking/tsmp_jube. Please ask for access if you want to use it.

With the *XML* file, the benchmark can be started by


1. Edit jube file for your needs
    * Load JUBE: ```module load JUBE```
    *    ```src_dir``` specify the absolute path to Cosmo, CLM, Oasis src files
    *    ```install_dir``` specify the path where to install TSMP
    *    provide/modify additional parameters if desired, e.g. account, testcase, etc.
2. Build TSMP executables
    * ```jube -vvv run ./TSMP_JUBE.xml --only-bench build```
3. To run the TSMP benchmark
    * ```jube -vvv run ./TSMP_JUBE.xml --only-bench run```
4. Make sure that all jobs are complete, e.g. with ```squeue -u $USER```
4. Proceed with analysis (modify ```bench_run``` if you specified a different name for the output directory)
    * Continue to execute jube script ```jube continue ./bench_run```
    * Analyse results ```jube analyse ./bench_run```
    * Show results ```jube result ./bench_run```
