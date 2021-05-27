# gain-variability

Implementation of gain variability estimation using a negative binomial distribution fit to spike count/rate, according to Goris et al., 2014 Nature Neuroscience doi:[10.1038/nn.3711](https://www.nature.com/articles/nn.3711).

Spiking activity is often modeled as a Poisson process in which spikes occur at random within a given time interval. One feature of the Poisson distribution is that it is parameterised by a single variable, representing both the variance and the mean of the distribution. 

However, the variance of spiking activity is often not equal to its mean, but larger. This super-Poisson variability can be modeled by fitting a negative binomial distribution instead of a Poisson distribution. The negative binomial distribution is parameterised by a dispersion parameter that captures this extra variability.

For an example on how to use this toolbox, see the [example application](https://htmlpreview.github.io/?https://github.com/jochemvankempen/gain-variability/blob/main/example.html).
