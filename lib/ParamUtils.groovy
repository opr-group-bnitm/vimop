// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


class ParamUtils {

    // CLI values are strings in Nextflow 26, so --x false would otherwise be truthy
    static boolean checkFlag(value) {
        return value != null && value.toString() != 'false'
    }

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
