# fms_distribution

Probability distributions

Every random variable, X, is determined by its _cumulative_distribution_function_. The function `cdf(x)` is the probability the random variable is less than `x`.

If a random variable takes on a finite or countable number of random variables it is _discrete_ and `pdf(x)` is the probability it takes on the value `x`. The function `is_discrete()` returns `true` and `domain()` returns the sequence of values it can take on.

A random variable is _continuously distributed_ if the cumulative distribution is absolutely continuous. The `pdf(x)` is the derivative of the cumulative distribution function at `x`. The function `is_continuous()` returns `true` and `domain()` returns a sequence of two values values: the lower and upper range of values at which the probablilty density is non-zero.

The _moment generating function_, E e<sup>t X</sup>, is given by `moment(t)`. The _moments_ are the sequence of Taylor coefficients of the power series and are returned by the function `moments()`

The _cumulant_, log E e<sup>s X</sup>, is given by `cumulant(t)`. The _cumulansts_ are the sequence of Taylor coefficients of the power series and are returned by the function `cumulants()`. Note the 0-th cumulant is always 0 and is not included in the sequence.
