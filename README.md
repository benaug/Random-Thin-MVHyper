# Random-Thin-MVHyper

A random thinning SCR model for subsampled detection data where you select X samples per trap-occasion. For example, for bear hair snares, it is common to select 1 sample per trap-occasion to genotype. The key assumption is that the sample selected is at random with respect to individual. This is a multivariate hypergeometric thinning model.

There are two observation models, each with and without a behavioral response to capture. The first is Poisson. This is a standard observation model for latent ID SCR, but it may not be realistic, particularly with a behavioral response to capture. In this latter case, individuals are both more likely to be captured and leave more samples conditional upon capture after a (positive) trap response. It seems more plausible that after a trap happy response, the capture probability goes up, but the number of samples left does not change. This is implied in the "hurdleZTPois" observation model. Here, there is a bernoulli detection process where the individual by trap detection probability is a function of the distance between the activity center and trap location. Then conditional on capture, individuals leave a zero-truncated Poisson number of samples. In the version with a behavioral response, only the detection probability undergoes a behavioral response.

See this paper for material on subsampling in the presence of a behavioral response to capture:

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12289

Then, see this paper where simulations show SCR seems robust to the effect identified above:

https://wildlife.onlinelibrary.wiley.com/doi/full/10.1002/jwmg.21144

But subsampling in this manner can introduce bias in the presence of a density gradient similar to here:

https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1748

If anyone wants to use this model and needs help, let me know--I think there is a paper to be written here.