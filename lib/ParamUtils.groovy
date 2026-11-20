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

    // Return a copy of `params` in which string values of integer, number and
    // boolean schema parameters are converted to those types, as Nextflow 24 did
    // for command-line values. Values that don't convert stay strings, so the
    // schema check still rejects them.
    static Map castToSchemaTypes(Map params, String schemaPath) {
        def schema = new groovy.json.JsonSlurper().parse(new File(schemaPath))
        Map types = [:]
        for (definition in schema['definitions'].values()) {
            for (property in definition['properties']) {
                types[property.key] = property.value['type']
            }
        }
        Map typed = new LinkedHashMap(params)
        for (name in params.keySet()) {
            typed[name] = castValue(params[name], types[name])
        }
        return typed
    }

    static Object castValue(value, type) {
        if (!(value instanceof CharSequence)) {
            return value
        }
        String text = value.toString().trim()
        if (type == 'boolean') {
            if (text.equalsIgnoreCase('true')) return true
            if (text.equalsIgnoreCase('false')) return false
        } else if (type == 'integer' || type == 'number') {
            if (text.isInteger()) return text.toInteger()
            if (text.isLong()) return text.toLong()
            if (type == 'number' && text ==~ /-?\d+(\.\d+)?([eE][-+]?\d+)?/) return text.toBigDecimal()
        }
        return value
    }
}
