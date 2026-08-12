/*
 * Shared parameter helpers.
 */

/*
 * Coerce a boolean-ish params value to a real boolean.
 *
 * Under the Nextflow ≥26.04 strict parser, CLI overrides like
 * `--filter_qc_failed false` reach the script as the String "false", which is
 * truthy in Groovy — so a plain `params.x ?` check silently reads it as "on".
 * 25.x coerced CLI booleans to Boolean before the script saw them. Comparing
 * the string form keeps one behavior across both parsers: exactly the values
 * true / "true" (any case) are on; everything else is off.
 */
def asBool(value) {
    return value.toString().toLowerCase() == 'true'
}
