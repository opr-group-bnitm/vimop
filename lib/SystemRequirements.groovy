// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.

import java.nio.file.FileSystems
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files
import java.lang.management.ManagementFactory
import com.sun.management.OperatingSystemMXBean
import nextflow.Nextflow


class SystemRequirements {
    boolean verbose

    SystemRequirements(boolean verbose = false) {
        this.verbose = verbose
    }

    void log(String message) {
        if (verbose) {
            println message
        }
    }

    static void exitError(String message) {
        def fullMessage = (
            "${message}\n"
            + "The minimum system requirements are an approximation and may not always fit.\n"
            + "You can change them by setting the respective parameters in EPI2ME (Miscellaneous Options), in the configs or via command line (type --help to list the options)."
        )
        System.err.println(fullMessage)
        System.exit(1)
    }

    static File findExistingParentDirectory(String dir) {
        def dirObj = Paths.get(dir).toAbsolutePath().toFile()
        while (dirObj != null && !dirObj.exists()) {
            dirObj = dirObj.getParentFile()
        }
        if(dirObj == null) {
            System.err.println("No existing parent directory found for $dir.")
            System.exit(1)
        }
        return dirObj
    }

    boolean checkDiskSpace(String dir, long minRequiredGB) {
        def file = findExistingParentDirectory(dir)
        def freeSpaceGB = file.usableSpace / (1024 * 1024 * 1024)
        log("Detected ${String.format('%.1f', freeSpaceGB)} GB free in directory ${dir} (Required: ${minRequiredGB} GB)")
        return freeSpaceGB >= minRequiredGB
    }

    boolean checkRam(long minRamGB) {
        try {
            OperatingSystemMXBean osBean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
            long totalRamGB = osBean.getTotalPhysicalMemorySize() / (1024 * 1024 * 1024)
            log("Detected system RAM: ${totalRamGB} GB (Required: ${minRamGB} GB)")
            return totalRamGB >= minRamGB
        } catch (Exception e) {
            System.err.println("Error checking RAM: ${e.message}")
            return false
        }
    }

    boolean checkCpus(int minCpus) {
        def availableCPUs = Runtime.runtime.availableProcessors()
        log("Detected ${availableCPUs} CPUs (Required: ${minCpus})")
        return availableCPUs >= minCpus
    }

    void checkSystemRequirements(
        long minDiskSpaceWorkGB,
        long minDiskSpaceOutGB,
        long minRamGB,
        int minCpus,
        String outputDir,
        String workDir
    ) {
        if (!checkDiskSpace(outputDir, minDiskSpaceOutGB)) {
            exitError("Insufficient disk space in output directory: $outputDir. Minimum required is ${minDiskSpaceOutGB}GB")
        }
        if (!checkDiskSpace(workDir, minDiskSpaceWorkGB)) {
            exitError("Insufficient disk space in work directory: $workDir. Minimum required is ${minDiskSpaceWorkGB}GB")
        }
        if (!checkRam(minRamGB)) {
            exitError("Insufficient RAM: Minimum required is ${minRamGB}GB")
        }
        if (!checkCpus(minCpus)) {
            exitError("Insufficient CPU cores: Minimum required is $minCpus")
        }
    }
}
