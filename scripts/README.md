# Miscellaneous script specs  

## `data_processing/make_variant_sets.py`  

Spec for `--sets-json` input is as follows:  

```
spec for .json input
{
  // Each set gets a unique identifier
  set_id:
  // Within each set is an array of sets of criteria.
  // The relationship between sets of criteria is assumed to be AND within 
  // each object (curly brackets) and assumed to be OR between objects
  [
    // Each key : value pair specifies a single criterion to require
    // The value is provided as an array of value and equality
    // Value can be numeric, boolean, or string
    // Equality should be two-letter string shorthand for comparison to apply
    // examples: eq, gt, lt, ge, le
    {
      key : [value, equality],
      key : [value, equality]
    },
    {
      key : [value, equality],
      key : [value, equality]
    }
  ]
}
```