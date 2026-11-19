// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


class ParamUtils {

    // Same rule as checkFlag() in main.nf. Values passed on the command line
    // arrive as strings under the strict syntax parser, so `--x false` would
    // otherwise be truthy.
    static boolean checkFlag(value) {
        return value != null && value.toString() != 'false'
    }
}
