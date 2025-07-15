Milani et al. Benchmark Data Generator
⚠️ Important Notice
Caution:
Execution of this tool will generate several hundred gigabytes of data. Ensure that the target storage location has sufficient free space. We strongly recommend at least 1.5 TB of available disk space.

Instructions
1. Prepare the Target Directory
Copy this entire repository folder to the (external) location where you intend to generate the output files.

2. Configure the Working Directory in Julia
Open the Julia script and ensure that the working directory is correctly set to the location where you have placed this folder.

Use the following commands in Julia to check or change the current directory:

julia
Copy
Edit
# Check the current working directory
pwd()

# Change the working directory
cd("desired/path/here")
Example:
If you want the data to be generated in:

bash
Copy
Edit
D:/datafiles/benchmark_data
Then the downloaded folder must be placed at:

bash
Copy
Edit
D:/datafiles/benchmark_data/Milani_et_al_benchmark_data
In Julia, check your current path:

julia
Copy
Edit
julia> pwd()
"C:\\Users\\analyst\\Documents\\Julia"
If this is not the intended location, change it:

julia
Copy
Edit
julia> cd("D:\\datafiles\\benchmark_data")
Verify the change:

julia
Copy
Edit
julia> pwd()
"D:\\datafiles\\benchmark_data"
Once the working directory is correctly set, proceed to the next step.

3. Run the Julia Script
Execute the Julia script in its entirety. Depending on your system's performance and the storage medium used, this process may take several hours to complete.

Recommendations:

Prevent your computer from entering sleep mode during execution.

Consider running the script overnight to avoid interruptions to other tasks.

4. Completion
After the script finishes running, all generated files will be located in your specified directory (e.g., D:/datafiles/benchmark_data in the example above).

You may then remove the Milani_et_al_benchmark_data folder from this location if it is no longer needed.

Notes
Ensure that Julia has the necessary permissions to read/write in the target directory.

Interrupting the script execution may result in incomplete data generation and require a rerun.

