<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
<h4 class="author">Mark Puttick</h4>
<h4 class="date">02 June 2019</h4>
</div>
<p>Models Of Trait Macroevolution On Trees (MOTMOT) is an R package that allows for testing of models of trait evolution (Thomas <em>et al.</em> 2012).</p>
<ul>
<li><a href="#models-of-trait-evolution">Tree transformation</a> models estimated using Maximum likelihood: <a href="#brownian-motion">Brownian motion</a>, <a href="#pagels-lambda">Pagel’s lambda</a>, <a href="#delta">Delta</a>, <a href="#kappa">Kappa</a>, <a href="#ornstein-uhlenbeck">Ornstein-Uhlenbeck (OU)</a>, <a href="#acdc-and-early-burst">Acceleration-Deaceleration (ACDC) and early bursts</a>, <a href="#psi-and-multispi">psi and multispi</a>, and <a href="#estimate-pagels-lambda-alongside-other-modes">estimating lambda alongside other models</a></li>
<li><a href="#rate-heterogeneous-models-of-evolution">Rate heterogeneous models of evolution</a>. Fit models in which the rate of evolution differs in clades selected <a href="#rate-heterogeneity-selected-a-priori"><em>a priori</em></a> (O’Meara <em>et al.</em> 2006; Thomas <em>et al.</em> 2006), and models with <a href="#rate-heterogeneity-with-no-a-priori-information">no <em>a-priori</em> shift locations</a> (Thomas <em>et al.</em> 2012)</li>
<li><a href="#timeslice-model">timeSlice</a> fit models in which all rates change at a specific time(s) by testing multiple shift times or those selected by the user</li>
<li><a href="#modeslice-model">modeSlice</a> fit models in which modes change at a specific time(s) in an extension to models introduced by Slater (2013)</li>
<li><a href="#nested-models-of-evolution">Nested Shift modes</a> Fit models models in which the ancestral BM rate switches to a ‘nested’ rate within a monophyletic clade in the phylogeny (Puttick 2018)</li>
<li><a href="#bayesian-estimation-of-tree-transformation-models">Bayesian estimation</a> of tree transformation models</li>
<li><a href="#character-displacement-models">Character displacement models</a> of inter-specific competition from Clarke <em>et al.</em> (2017)</li>
<li><a href="#fast-estimation-of-phylogenetic-generalised-least-squares">Fast estimation of Phylogenetic Generalised Least Squares (PGLS)</a> using independent contrasts</li>
</ul>
<div id="introduction" class="section level1">
<h1>Introduction</h1>
<p>First we install</p>
<pre class="r"><code>install.packages(&quot;motmot&quot;)</code></pre>
<p>and load MOTMOT</p>
<pre class="r"><code>library(motmot)</code></pre>
<p>For these examples we will use anoles lizard data available from MOTMOT. A time-calibrated phylogeny of <em>Anolis</em> species <code>anolis.tree</code>, and various trait and biogeographical trait data <code>anolis.data</code>.</p>
<pre class="r"><code>data(anolis.tree)
data(anolis.data)
attach(anolis.data)
anolis.tree</code></pre>
<pre><code>## 
## Phylogenetic tree with 165 tips and 164 internal nodes.
## 
## Tip labels:
##  A_occultus, A_darlingt, A_monticol, A_bahoruco, A_dolichoc, A_henderso, ...
## Node labels:
##  2, 2, 2, 2, 2, 2, ...
## 
## Rooted; includes branch lengths.</code></pre>
<p>We will use the continuous trait data: male snout-ventral length <code>Male_SVL</code>. Here, we construct a matrix of just <code>Male_SVL</code> data, remove missing data, and log-transform the values. All this can be done using the function <code>sortTraitData</code></p>
<pre class="r"><code>sortedData &lt;- sortTraitData(phy = anolis.tree, y = anolis.data, 
    data.name = &quot;Male_SVL&quot;, pass.ultrametric = TRUE)
phy &lt;- sortedData$phy
male.length &lt;- sortedData$trait</code></pre>
<p>Finally, we will ‘prune’ the species from the tree using <code>drop.tip</code> from <a href="https://CRAN.R-project.org/package=ape">APE</a>. We plot our tree and data using the MOTMOT <code>traitData.plot</code> function.</p>
<pre class="r"><code>traitData.plot(y = male.length, phy, lwd.traits = 2, col.label = &quot;#00008050&quot;, 
    tck = -0.01, mgp = c(0, 0.2, 0), cex.axis = 0.5, show.tips = FALSE)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot1-1.png" alt="Figure 1. TraitData showing the realtive male snout-vent length at the tips" width="1000" />
<p class="caption">
Figure 1. TraitData showing the realtive male snout-vent length at the tips
</p>
</div>
<p>For the sake of brevity, in the following examples we fit the models to a subset of these data: including the clade from node 182 only using the <a href="https://CRAN.R-project.org/package=ape">APE</a> function <code>extract.clade</code>.</p>
<pre class="r"><code>## uncomment to view the tree plot(phy, show.tip.label=FALSE,
## no.margin=TRUE, edge.col=&#39;grey20&#39;) nodelabels(182, 182,
## bg=&#39;black&#39;, col=&#39;white&#39;)
phy.clade &lt;- extract.clade(phy, 182)
male.length.clade &lt;- as.matrix(male.length[match(phy.clade$tip.label, 
    rownames(male.length)), ])</code></pre>
</div>
<div id="models-of-trait-evolution" class="section level1">
<h1>Models of trait evolution</h1>
<p>We can now test various models of evolution using our trait data.</p>
<div id="brownian-motion" class="section level2">
<h2>Brownian motion</h2>
<p>To start we will fit a simple Brownian motion model to the data, as the null hypothesis of phylogenetic trait evolution (Cavalli-Sforza and Edwards 1967; Felsenstein 1973; 1985). Brownian motion describes a process in which tip states are modelled under the assumption of a multi-variate normal distribution. On a phylogeny, the multi-variate mean of tip states is equal to the root state estimate, and variance accummulates linearly through time. Trait evolution is shared but following a split individual branches evolve and accummulate trait variance independently from their shared ancestral value.</p>
<p>The function <code>transformPhylo.ML</code> is used to fit Brownian motion models and its derivatives. Here we fit a simple Brownian motion model to the subset of anolis male SVL data to obtain the Brownian variance, ancestral estimate, log-likelihood, Akaike Information Criterion (AIC), and small-sample AIC (AICc).</p>
<pre class="r"><code>bm.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;bm&quot;)
bm.ml</code></pre>
<pre><code>## $brownianVariance
##             [,1]
## [1,] 0.001858067
## 
## $logLikelihood
## [1] -0.2838382
## 
## $root.state
## [1] 3.849481
## 
## $AIC
## [1] 4.567676
## 
## $AICc
## [1] 5.047676</code></pre>
</div>
<div id="pagels-lambda" class="section level2">
<h2>Pagel’s lambda</h2>
<p>Here we fit models to test Pagel’s lambda (Pagel 1997; 1999). Pagel’s lambda is a measure of phylogenetic ‘signal’ in which the degree to which shared history of taxa has driven trait distributions at tips. In this model, internal branch lengths are transformed by the lambda parameter value. When the parameter lambda equals 1, branches are transformed by multiplying by 1 and so the model is equal to Brownian motion (high phylogenetic signal). Values of lambda under 1 suggest there has been less influence of shared history on trait values at the tips. Finally, a lambda value of 0 indicates no phylogenetic influence on trait distributions, and is equivalent to a ‘star phylogeny’ with no shared branch lengths.</p>
<p>The maximum likelihood of lambda can be estimated in MOTMOT.</p>
<pre class="r"><code>lambda.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;lambda&quot;)
lambda.ml</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 6.522191
## 
## $Lambda
##       MLLambda   LowerCI   UpperCI
## [1,] 0.8369973 0.5259423 0.9742338
## 
## $brownianVariance
##              [,1]
## [1,] 0.0008245357
## 
## $root.state
## [1] 3.853432
## 
## $AIC
## [1] -7.044383
## 
## $AICc
## [1] -6.044383</code></pre>
<p>The maximum likelhood estimate of Pagel’s lambda is equal to 0.84.</p>
<p>A new feature in MOTMOT allows for plotting of the likelihood profile for the branch-transformation parameter, in this case Pagel’s lambda using the argument <code>profilePlot</code> in <code>transformPhylo.ML</code>.</p>
<pre class="r"><code>lambda.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;lambda&quot;, profilePlot = TRUE)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot2-1.png" alt="Figure 2. Profile plot of ML estimation for Pagel&#39;s lambda" width="1000" />
<p class="caption">
Figure 2. Profile plot of ML estimation for Pagel’s lambda
</p>
</div>
<p>We now compare the relative fit of the BM and lambda models. Lambda has higher likelihood, but it also has more parameters. The root state and sigma-squared (rate) parameters are present in both models but the lambda model also requires an estimate of the parameter lambda. We can test whether the lambda model is a significant improvement over BM. First we test the relative fit by using the chi-squared distribution. The models differ in one degree of freedom: BM has 2 parameters and lambda has 3. We can use the <code>stats</code> function <code>pchisq</code> to obtain a p value by testing using a chi-squared distribution. The lambda is indeed a superior fit compared to BM when fit to these data (<em>p</em> &lt; 0.05).</p>
<pre class="r"><code>p.value &lt;- 1 - pchisq(lambda.ml$MaximumLikelihood - bm.ml$logLikelihood, 
    1)
p.value</code></pre>
<pre><code>## [1] 0.009085056</code></pre>
<p>Additionally there is a large small-sample Akaike Information Criterion (AICc) difference between the two models: BM has a higher AICc compared to lambda. The difference (11.09) is &gt;4 which is traditionally seen as indication of a superior fit (Burnham and Anderson 2003).</p>
<pre class="r"><code>bm.ml$AICc - lambda.ml$AICc</code></pre>
<pre><code>## [1] 11.09206</code></pre>
</div>
<div id="delta" class="section level2">
<h2>Delta</h2>
<p>Delta indicates a decrease or increase in the rate of trait evolution through time (Pagel 1997; 1999); a value of 1 is equivalent to Brownian motion, &lt; 1 indicates a slow-down, and &gt; 1 is indicates greater change closer to the present. Here we find a Maximum likelihood estimated for Delta of 2.23 but the CI spans &lt; 1 to &gt; 4, so it is not possible to conclusively support a change in the rate of evolution through time.</p>
<pre class="r"><code>delta.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;delta&quot;)
delta.ml</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 1.179797
## 
## $Delta
##       MLDelta   LowerCI  UpperCI
## [1,] 2.231464 0.8477531 4.660712
## 
## $brownianVariance
##              [,1]
## [1,] 6.013993e-06
## 
## $root.state
## [1] 3.8843
## 
## $AIC
## [1] 3.640407
## 
## $AICc
## [1] 4.640407</code></pre>
</div>
<div id="kappa" class="section level2">
<h2>Kappa</h2>
<p>Kappa is used as a measure of punctuated evolution and spans values of 0-1 (Pagel 1997:1999). A Kappa value of 1 is equivalent to BM, and 0 indicates trait change occurs at events of speciation. Here there is evidence of punctuated evolution. <code>transformPhylo.ML</code> also allows users to see the the phylogeny transformed by model parameters. As an example, we show the original, BM model phylogeny and compare this with the phylogeny transformed by the Kappa phylogeny.</p>
<pre class="r"><code>kappa.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;kappa&quot;, profilePlot = FALSE, returnPhy = TRUE)
par(mfrow = c(1, 2))
plot.phylo(phy.clade, show.tip.label = FALSE, no.margin = TRUE)
mtext(&quot;Original phylogeny&quot;, 3, cex = 0.7, line = -1)
plot.phylo(kappa.ml$kappaPhy, show.tip.label = FALSE, no.margin = TRUE)
mtext(&quot;Kappa model phylogeny&quot;, 3, cex = 0.7, line = -1)
mtext(&quot;Kappa = 1e-8&quot;, 3, cex = 0.7, line = -2)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot3-1.png" alt="Figure 3. Comparison of BM and Kappa transformed trees." width="1000" />
<p class="caption">
Figure 3. Comparison of BM and Kappa transformed trees.
</p>
</div>
</div>
<div id="ornstein-uhlenbeck" class="section level2">
<h2>Ornstein-Uhlenbeck</h2>
<p>The OU model allows for modelling of attraction to a optimum value, alpha (Hansen 1997; Butler and King 2004). This model again is similar to the Brownian motion model, but models the strength of attraction to alpha. The OU model can be difficult to interpret and care is advised in its use (Cooper <em>et al.</em> 2016).</p>
<p>In MOTMOT, as with most implements of the phylogenetic OU model, the value of the optimum is equal to the ancestral trait estimate. With all <code>transformPhylo.ML</code> functions it is possible to change the bounds on the estimated parameters. For example, here the value of <em>alpha</em> is constrained to 2 using the argument <code>upperBound</code>.</p>
<pre class="r"><code>ou.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;OU&quot;, profilePlot = TRUE, upperBound = 2)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot4-1.png" alt="Figure 4. Profile plot to estimate alpha" width="1000" />
<p class="caption">
Figure 4. Profile plot to estimate alpha
</p>
</div>
<pre class="r"><code>ou.ml</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 1.743459
## 
## $Alpha
##         MLAlpha      LowerCI    UpperCI
## [1,] 0.01693353 0.0004499003 0.03912232
## 
## $brownianVariance
##             [,1]
## [1,] 0.002823932
## 
## $root.state
## [1] 3.876702
## 
## $AIC
## [1] 2.513082
## 
## $AICc
## [1] 3.513082</code></pre>
<p>The value of alpha is higher than zero, but very small (0.01692855). So the model is not equivalent to Brownian motion but there is little evidence from AICc that the model is an improvement, and the likelihood ratio test show a non-significant improvement so it does not have higher relative support compared to BM (<em>p</em> &gt; 0.05).</p>
<pre class="r"><code>p.value &lt;- 1 - pchisq(ou.ml$MaximumLikelihood - bm.ml$logLikelihood, 
    1)
p.value</code></pre>
<pre><code>## [1] 0.1544952</code></pre>
<pre class="r"><code>bm.ml$AICc - ou.ml$AICc</code></pre>
<pre><code>## [1] 1.534594</code></pre>
</div>
<div id="acdc-and-early-burst" class="section level2">
<h2>ACDC and Early Burst</h2>
<p>A new addition to MOTMOT is the ACDC model (Blomberg <em>et al.</em> 2003). This model allows for exponential changes in the rate of evolution in the history of a clade.</p>
<pre class="r"><code>acdc.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;ACDC&quot;, profilePlot = TRUE)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot5-1.png" alt="Figure 5. Profile plot to estimate the ACDC parameter" width="1000" />
<p class="caption">
Figure 5. Profile plot to estimate the ACDC parameter
</p>
</div>
<pre class="r"><code>acdc.ml</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 1.743459
## 
## $ACDC
##          MLacdc     LowerCI    UpperCI
## [1,] 0.03385981 0.000879245 0.04246516
## 
## $brownianVariance
##              [,1]
## [1,] 0.0002590746
## 
## $root.state
## [1] 3.876696
## 
## $AIC
## [1] 2.513082
## 
## $AICc
## [1] 3.513082</code></pre>
<p>There is little evidence here of exponential decreases or increases in the rate of trait evolution - the ACDC exponential parameter is close to 0 (0.034). We can see this is not a significant improvement on BM.</p>
<pre class="r"><code>p.value.2 &lt;- 1 - pchisq(acdc.ml$MaximumLikelihood - bm.ml$logLikelihood, 
    1)
p.value.2</code></pre>
<pre><code>## [1] 0.1544951</code></pre>
<p>As an example, here we constrain the ‘upperBound’ to &lt; 0, this is equivalent to the Early Burst model (Harmon <em>et al.</em> 2010) fit in <a href="https://CRAN.R-project.org/package=geiger">geiger</a>.</p>
<pre class="r"><code>transformPhylo.ML(phy = phy.clade, y = male.length.clade, model = &quot;ACDC&quot;, 
    profilePlot = FALSE, upperBound = -1e-06, print.warning = FALSE)</code></pre>
<pre><code>## $MaximumLikelihood
## [1] -0.2839606
## 
## $ACDC
##      MLacdc     LowerCI UpperCI
## [1,] -1e-06 -0.01322547  -1e-06
## 
## $brownianVariance
##             [,1]
## [1,] 0.001858181
## 
## $root.state
## [1] 3.849481
## 
## $AIC
## [1] 6.567921
## 
## $AICc
## [1] 7.567921</code></pre>
<p>The estimate of -1e-6 for the exponential decrease parameter, which means the model is effectively equivalent to Brownian motion.</p>
</div>
<div id="psi-and-multispi" class="section level2">
<h2>psi and multispi</h2>
<p>The parameter psi is similar to the parameter Kappa in that it is a measure of the relative contribution of speciational (~punctuated) and gradual evolution to trait change on a phylogeny (Ingram 2011; Ingram <em>et al.</em> 2016). The parameter psi is based upon measures of evolution over time and at speciation, and can also account for ‘hidden’ nodes not seen in the input phylogeny. The parameter psi measures the proportion of total evolutionary change (speciational + gradual) that can be attributable to speciational evolution, so the estimation for psi between 0 (Brownian motion) and 1 (indicating equal branch lengths, ~speciational change).</p>
<p>In MOTMOT we can fit a simple psi model using the input tree.</p>
<pre class="r"><code>psi.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;psi&quot;, profilePlot = TRUE)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot6-1.png" alt="Figure 6. Profile plot to estimate the psi parameter" width="1000" />
<p class="caption">
Figure 6. Profile plot to estimate the psi parameter
</p>
</div>
<pre class="r"><code>psi.ml</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 8.495569
## 
## $psi
##      MLpsi  LowerCI UpperCI
## [1,]     1 0.316058       1
## 
## $brownianVariance
##              [,1]
## [1,] 0.0008018518
## 
## $root.state
## [1] 3.761269
## 
## $AIC
## [1] -10.99114
## 
## $AICc
## [1] -9.991139</code></pre>
<p>This indicates support for the psi model is a significant improvement on Brownian motion (<em>p</em> &lt; 0.05).</p>
<pre class="r"><code>p.value.psi &lt;- 1 - pchisq(psi.ml$MaximumLikelihood - bm.ml$logLikelihood, 
    1)
p.value.psi</code></pre>
<pre><code>## [1] 0.003046501</code></pre>
<p>We could also get a potentially more accurate of speciation rates by using the full tree, rather than the pruned tree to estimate speication and extinction rates as this will give more accurate estimates rather than using the taxa with complete data only. If extinction rates are larger than 0, then the estimates will differ from the simple model above.</p>
<pre class="r"><code>psi_ext.est &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;psi&quot;, profilePlot = FALSE, hiddenSpeciation = TRUE, 
    full.phy = phy)
all.equal(psi.ml, psi_ext.est)</code></pre>
<pre><code>## [1] TRUE</code></pre>
<p>In this case, there is no difference in the estimates as extinction rates are equal to 0.</p>
<p>We can also apply multipsi model in which different regions of the tree have different estimates of the parameter psi. We can now fit the multispi model with these data. In MOTMOT, this model requires branch labels given <em>a priori</em> by the user to delimit the different regimes on the phylogeny. Note that these clades with potentially different psi regimes do not need to be monophyletic clades. Here we arbitarily assign two clades ‘a’ and ‘b’ to test differences between them.</p>
<pre class="r"><code>plot(phy.clade, no.margin = TRUE, cex = 0.8)
two.clade.labels &lt;- c(rep(&quot;a&quot;, 17), rep(&quot;b&quot;, 37))
edgelabels(two.clade.labels, col = c(rep(&quot;blue&quot;, 17), rep(&quot;red&quot;, 
    37)), bg = &quot;white&quot;)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot7-1.png" alt="Figure 7. Two clades used in the multipsi model" width="1000" />
<p class="caption">
Figure 7. Two clades used in the multipsi model
</p>
</div>
<p>Using these data we fit the model with <code>transformPhylo.ML</code>.</p>
<pre class="r"><code>transformPhylo.ML(phy = phy.clade, y = male.length.clade, model = &quot;multipsi&quot;, 
    branchLabels = c(rep(&quot;a&quot;, 17), rep(&quot;b&quot;, 37)), hiddenSpeciation = TRUE, 
    full.phy = phy)</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 8.495569
## 
## $psi
##   MLpsi    LowerCI UpperCI
## a     1 0.04812257       1
## b     1 0.20955316       1
## 
## $brownianVariance
##              [,1]
## [1,] 0.0008018518
## 
## $root.state
## [1] 3.761269
## 
## $AIC
## [1] -8.991139
## 
## $AICc
## [1] -7.252008</code></pre>
<p>In this model, the estimate of psi does not differ between the two regions of the tree</p>
</div>
<div id="estimate-pagels-lambda-alongside-other-modes" class="section level2">
<h2>Estimate Pagel’s lambda alongside other modes</h2>
<p>One way to deal with ‘noisy’ data is to estimate Pagel’s lambda alongside a parameter of interest. By using Pagel’s lambda alongside other models it may be possible to account for variation in the data that may be a result of errors in the phylogeny or trait data. In MOTMOT, Pagel’s lambda can be estimated alongside the delta, kappa, OU, psi, and ACDC models. Here we look at example using ACDC. The model is fit with same function. <code>transformPhyo.ML</code> but with the argument <code>lambdaEst</code> set to <code>TRUE</code>.</p>
<pre class="r"><code>acdc.ml.lambda &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;ACDC&quot;, lambdaEst = TRUE)
# original ACDC model
acdc.ml</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 1.743459
## 
## $ACDC
##          MLacdc     LowerCI    UpperCI
## [1,] 0.03385981 0.000879245 0.04246516
## 
## $brownianVariance
##              [,1]
## [1,] 0.0002590746
## 
## $root.state
## [1] 3.876696
## 
## $AIC
## [1] 2.513082
## 
## $AICc
## [1] 3.513082</code></pre>
<pre class="r"><code># ACDC model plus lambda
acdc.ml.lambda</code></pre>
<pre><code>## $MaximumLikelihood
## [1] 7.376867
## 
## $ACDC
##          MLacdc    LowerCI     UpperCI
## [1,] -0.1829847 -0.3244672 -0.08291527
## 
## $brownianVariance
##            [,1]
## [1,] 0.01272729
## 
## $root.state
## [1] 3.83604
## 
## $lambda
## [1] 0.1388712
## 
## $AIC
## [1] -4.753735
## 
## $AICc
## [1] -2.026462</code></pre>
<p>We can see lambda is &lt; 1, and this has affected the parameter estimation. The improvement in the model fit is significant compared to the ACDC model fit without estimating lambda and the null BM model.</p>
<pre class="r"><code># p value of the ACDC and ACDC+lambda models. No significant
# improvement
1 - pchisq(acdc.ml.lambda$MaximumLikelihood - acdc.ml$MaximumLikelihood, 
    df = 1)</code></pre>
<pre><code>## [1] 0.01762134</code></pre>
<pre class="r"><code># p value of the BM and ACDC+lambda model comparison. No
# significant improvement
1 - pchisq(acdc.ml.lambda$MaximumLikelihood - bm.ml$logLikelihood, 
    df = 2)</code></pre>
<pre><code>## [1] 0.02170196</code></pre>
</div>
</div>
<div id="rate-heterogeneous-models-of-evolution" class="section level1">
<h1>Rate heterogeneous models of evolution</h1>
<div id="rate-heterogeneity-selected-a-priori" class="section level2">
<h2>rate heterogeneity selected <em>a priori</em></h2>
<p>MOTMOT can test models of evolution in which pre-defined clades can vary in the rate of evolution. Here we fit a model in which the nodes descending from nodes 32 and 49 have a seperate rate of evolution. First, we can visualise these nodes on the phylogeny.</p>
<pre class="r"><code>plot(phy.clade, show.tip.label = FALSE, no.margin = TRUE, edge.col = &quot;grey20&quot;)
nodelabels(c(32, 49), c(32, 49), bg = &quot;black&quot;, col = &quot;white&quot;)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot8-1.png" alt="Figure 8. Lineages with different rates of evolution" width="1000" />
<p class="caption">
Figure 8. Lineages with different rates of evolution
</p>
</div>
<p>We then fit the MOTMOT model, again using the function <code>transformPhylo.ML</code>. We use the argument <code>model=clade</code>. This fits the non-censored model of O’Meara <em>et al.</em> (2006).</p>
<pre class="r"><code>cladeRate.ml &lt;- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
    model = &quot;clade&quot;, nodeIDs = c(32, 49))
cladeRate.ml</code></pre>
<pre><code>## $Rates
##      node    MLRate   LowerCI  UpperCI
## [1,]   32 0.8138100 0.2632955 3.263628
## [2,]   49 0.6819079 0.1896347 2.952364
## 
## $MaximumLikelihood
## [1] -0.1462557
## 
## $brownianVariance
##             [,1]
## [1,] 0.002143258
## 
## $root.state
## [1] 3.870488
## 
## $AIC
## [1] 8.292511
## 
## $AICc
## [1] 10.03164</code></pre>
<p>These results indicate that the two clades tend to have a lower rate of evolution compared to the background rate. However, the CIs indicate these decreases may not be robust.</p>
</div>
<div id="rate-heterogeneity-with-no-a-priori-information" class="section level2">
<h2>rate heterogeneity with no <em>a priori</em> information</h2>
<p>We can also fit rate heterogeneous models without specifying where we expect shifts on the tree. We can use the arguments <code>model=&quot;tm1&quot;</code> and <code>model=&quot;tm2&quot;</code>. These models fit the <code>traitMedusa</code> model in which nodes are individually tested for rate increases or decreases (Thomas and Freckleton 2012), and the clade or branch with a rate change that produces the largest increase in likelihood is returned. Note, it is possible to exclude small nodes using the argument <code>minCladeSize</code>. As well as allowing clade differences in rate, the <code>tm2</code> also allows for branch-based increases or decreases in rate, whereas <code>tm1</code> only searches for clade-based rate changes.</p>
<p>We can now fit the <code>tm2</code> algorithm. The output shows the log-likelihood, AIC, AICc, rate type (branch of clade), for the best-fitting model at each stage. This starts with the BM model, and then one shift model, two shift model, etc.,</p>
<pre class="r"><code># tm1 algorithm not run tm1.ml &lt;-
# transformPhylo.ML(y=male.length.clade, phy=phy.clade,
# model=&#39;tm1&#39;, minCladeSize=2, nSplits=3)
# trait.medusa.tm1.summary &lt;- traitMedusaSummary(tm1.ml,
# cutoff=2, AICc=T) tm2 model
tm2.ml &lt;- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
    model = &quot;tm2&quot;, minCladeSize = 5, nSplits = 2)</code></pre>
<pre><code>## 
##  BM model
##     node shiftPos        lnL n.params      AIC     AICc
## BM     0        1 -0.2838382        2 4.567676 5.047676
## 
##  Shift 1
##         node shiftPos      lnL n.params         AIC     AICc     rate.1
## shift.1   39    clade 3.042358        3 -0.08471602 0.915284 0.09148646
## 
##  Shift 2
##         node shiftPos      lnL n.params       AIC     AICc    rate.1
## shift.2   44    clade 4.746785        5 0.5064296 3.233702 0.1408068
##           rate.2
## shift.2 3.158565</code></pre>
<p>We can summarise the analyses using <code>traitMedusaSummary</code> and plotting the shifts on the phylogeny using the function <code>plotPhylo.motmot</code>. These results show a decrease at node 39 that we can visualise on the phylogeny.</p>
<pre class="r"><code>trait.medusa.tm2.summary &lt;- traitMedusaSummary(tm2.ml, cutoff = 2, 
    AICc = TRUE)
trait.medusa.tm2.summary</code></pre>
<pre><code>## $ModelFit
##              lnL n.params         AIC     AICc
## shift.1 3.042358        3 -0.08471602 0.915284
## 
## $Rates
##   node shiftPos             MLRate            LowerCI           UpperCI
## 1   39    clade 0.0914864604723702 0.0260301808222033 0.500374159059036
## 
## $optimalTree
## 
## Phylogenetic tree with 28 tips and 27 internal nodes.
## 
## Tip labels:
##  A_alutaceu, A_inexpect, A_vanidicu, A_alfaroi, A_macilent, A_clivicol, ...
## Node labels:
##  2, 2, 2, 2, 2, 2, ...
## 
## Rooted; includes branch lengths.</code></pre>
<pre class="r"><code>colour_motmot &lt;- plotPhylo.motmot(phy = phy.clade, traitMedusaObject = trait.medusa.tm2.summary, 
    reconType = &quot;rates&quot;, type = &quot;fan&quot;, cex = 0.5, edge.width = 2)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot9-1.png" alt="Figure 9. The subset of the tree showing the rate heterogeneity estimated from the traitMedusa model" width="1000" />
<p class="caption">
Figure 9. The subset of the tree showing the rate heterogeneity estimated from the traitMedusa model
</p>
</div>
<p>Thomas and Freckleton (2012) showed the <code>tm2</code> algortihm has a high type-one error rate. One way to ameliorate this is to estimate the level a one shift is supported when we know BM is the true model. For example, we could simulate 1000 BM datasets on the tree, estimate a single shift using the <code>tm2</code> algortihm, and calculating the difference between the AICcs for each BM and one shift model. We can these use this difference to estimate the AICc ‘penalty’ the is needed to reduce the <code>tm2</code> type-one error rate to 0.05. We could use this penalty in the <code>cutoff</code> argument of the <code>traitMedusaSummary</code> argument.</p>
<p>This can all be calculated with the MOTMOT function <code>calcCutOff</code>. The function requires the tree and input from the model applied to the empirical data as well as the number of simulations. Here we calculated the AICc cut-off required for the <code>tm2</code> analysis from above (for brevity this is not run here, but should be run for each analysis individually).</p>
<pre class="r"><code>## uncomment to run set.seed(203); calcCutOff(phy.clade,
## n=1000, model=&#39;tm2&#39;, minCladeSize=5, nSplits=1); 95%
## 5.698198</code></pre>
<p>Here if we repeat this analysis with the appropriate AICc cut-off (5.698) the we see that the single-rate Brownian motion is, in fact, supported.</p>
<pre class="r"><code>traitMedusaSummary(tm2.ml, cutoff = 5.698198, AICc = TRUE)$Rates</code></pre>
<pre><code>## [1] &quot;Single rate&quot;</code></pre>
</div>
</div>
<div id="timeslice-model" class="section level1">
<h1>timeSlice model</h1>
<p>A new addition to motmot is a Maximum likelihood model that allows for heterogeneous rates in different time periods. These models are seperate from the models that allow for heterogeneous rates among lineages, as modelled by the <code>traitMedusa</code> algorithms.</p>
<p>The <code>timeSlice</code> model is implemented using the <code>transformPhylo.ML</code> function, using the argument <code>model=&#39;timeSlice&#39;</code>. The function allows for two seperate models of evolution. In one, it is possible to test shifts in evolution at times selected <em>a priori</em>. Alternatively, the fit of models can be tested at a range of different times, and the function will return the best-fitting model</p>
<p>First we will test for a shift in the rate of evolution 10 million years ago.</p>
<pre class="r"><code>timeSlice.10.ml &lt;- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
    model = &quot;timeSlice&quot;, splitTime = 10)</code></pre>
<pre><code>## [1] &quot;BM model&quot;
##          lnL          AIC         AICc   sigma.sq.1  anc.state.1 
## -0.283838169  4.567676338  5.047676338  0.001858067  3.849481405 
## [1] &quot;shiftModel&quot;
##          lnL          AIC         AICc   sigma.sq.1  anc.state.1 
##  2.946497537  2.107004926  3.846135361  0.001006387  3.860015270 
##       rates1       rates2  time.split1 
##  0.692073080  2.944764076 10.000000000</code></pre>
<p>We can use the function <code>timeSliceSummary</code> to summarise and plot the results. The output summarises the best model according to AICc fit. This function automatically plots the original tree showing the location of shift(s), and the colours show the relative rates in each time slice. The second plot below shows the same tree and colours, but with the branch lengths scaled to the ML optimised rates</p>
<pre class="r"><code>outputSummary &lt;- timeSliceSummary(timeSlice.10.ml, cutoff = 0.001, 
    cex = 0.55, edge.width = 2, cex.plot = 0.8, colour.ramp = c(&quot;blue&quot;, 
        &quot;red&quot;), label.offset = 0.5)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot10-1.png" alt="Figure 10. TimeSlice plot with a split at 10 Ma" width="1000" />
<p class="caption">
Figure 10. TimeSlice plot with a split at 10 Ma
</p>
</div>
<p>We can also see other summarise information, such as the CI for each rate estimate.</p>
<pre class="r"><code>outputSummary$RatesCI</code></pre>
<pre><code>##            rates       LCI      UCI
## rates1 0.6920731 0.1873339  2.11563
## rates2 2.9447641 0.9633084 10.87739</code></pre>
<p>Rather than testing the overall fit of each model, the model can search all shift times and returns the shift location or locations with the highest likelihood. The function automatically tests for all 1 Ma shifts between the age of the tree - 10 Ma, and the present + 10 Ma; all these presets can be customised using the <code>boundaryAge</code> argument that supplies a vector with the first age specifying the distance from the root and the second age specifying the age from the tips. The <code>splitTime</code> argument sets the ages at which all shifts will be tested for between the <code>boundaryAge</code> with the default testing all shifts at 1 Ma intervals. The model searches for <em>n</em> shifts set by the <code>nSplits</code> argument.</p>
<p>This model searches for the highest likelihood single shift by searching for the highest likelihood shift time between 62-8 Myrs.</p>
<pre class="r"><code>timeSlice.ml &lt;- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
    model = &quot;timeSlice&quot;, nSplits = 1, boundaryAge = 8)</code></pre>
<pre><code>## [1] BM model
##          lnL          AIC         AICc   sigma.sq.1  anc.state.1 
## -0.283838169  4.567676338  5.047676338  0.001858067  3.849481405 
## [1] shift 1
##         lnL         AIC        AICc  sigma.sq.1 anc.state.1      rates1 
## 3.584675298 0.830649404 2.569779838 0.001012562 3.859446535 0.666861158 
##      rates2 time.split1 
## 3.238675325 8.000000000</code></pre>
<p>And summarise the results. We can selected the cutoff AICc improvement needed to justify selecting the next model. Here we use the arbitary cut-off value of 1. We could test this formally by estimating the correct AICc value needed to reduced type-error &gt; 5% by using BM simulated data (an example using the tm2 is shown above)</p>
<pre class="r"><code>outputSummary &lt;- timeSliceSummary(timeSlice.ml, cutoff = 1, cex = 0.2, 
    edge.width = 2, cex.plot = 0.8, colour.ramp = c(&quot;blue&quot;, &quot;red&quot;), 
    label.offset = 0.5)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot11-1.png" alt="Figure 11. TimeSlice plot with Maximum likelihood estimation of split time" width="1000" />
<p class="caption">
Figure 11. TimeSlice plot with Maximum likelihood estimation of split time
</p>
</div>
</div>
<div id="modeslice-model" class="section level1">
<h1>modeSlice model</h1>
<p>In a related extension, we have incorporated the new <code>modeSlice</code> model to the <code>transformPhylo.ML</code>. <code>modeSlice</code> incorporates and extends the methods of Slater (2013) by allowing for multiple shifts in various modes of evolution (BM, OU, EB, and Kappa) at different times in the phylogeny’s history. This is flexible as users can input multiple rate shift times with different combinations of modes. Furthermore, time bins with a BM mode of evolution can optionally vary in the rate of evolution compared to the background variance (<code>rate.var</code> argument), and users can include a rate scalar alongside EB modes.</p>
<p>Here a model is fit with a shift from an EB model with associated rate scalar to an OU model 40 Ma and then to a BM rate shift model at 30 Ma to the present. The results indicate an ACDC/EB scalar (root age-40Ma), followed by a OU model with alpha of 1.75 (40-30Ma), followed by a rate increase (5.3x background from 30-0 Ma). However this model is not supported over Brownian motion.</p>
<pre class="r"><code>modeSlice.ml &lt;- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
    model = &quot;modeSlice&quot;, splitTime = c(40, 30), mode.order = c(&quot;ACDC&quot;, 
        &quot;OU&quot;, &quot;BM&quot;), rate.var = TRUE, acdcScalar = TRUE)
modeSlice.ml$AICc</code></pre>
<pre><code>## [1] 14.72642</code></pre>
<pre class="r"><code>bm.ml$AICc</code></pre>
<pre><code>## [1] 5.047676</code></pre>
</div>
<div id="nested-models-of-evolution" class="section level1">
<h1>Nested models of evolution</h1>
<p>We can also tested models of nested evolution in which an ancestral model of BM evolution changes to a alternative model (EB, OU, kappa, delta, psi) within the phylogeny (Puttick 2018).</p>
<p>Here we can show an example of BM -&gt; OU and BM -&gt; ACDC at node 44 of the phylogeny. However, neither of these is a significantly better relative fit than BM.</p>
<pre class="r"><code>bm.model &lt;- transformPhylo.ML(male.length.clade, phy = phy.clade, 
    model = &quot;bm&quot;)
nested.acdc &lt;- transformPhylo.ML(male.length.clade, phy = phy.clade, 
    model = &quot;ACDC&quot;, nodeIDs = 44)
nested.ou &lt;- transformPhylo.ML(male.length.clade, phy = phy.clade, 
    model = &quot;OU&quot;, nodeIDs = 44)
1 - pchisq(nested.acdc$MaximumLikelihood - bm.model$logLikelihood, 
    1)</code></pre>
<pre><code>## [1] 0.05740889</code></pre>
<pre class="r"><code>1 - pchisq(nested.ou$MaximumLikelihood - bm.model$logLikelihood, 
    1)</code></pre>
<pre><code>## [1] 0.3614244</code></pre>
</div>
<div id="bayesian-estimation-of-tree-transformation-models" class="section level1">
<h1>Bayesian estimation of tree transformation models</h1>
<p>Parameters of various modes of evolution can be conducted using a simple Bayesian Markov Chain Monte Carlo (MCMC) algorithm in <code>transformPhylo.MCMC</code> which may better reflect probabilistic uncertainty in parameter estimates compared to Maximum likelihood estimation. By default the model places a uniform prior and new proposals using an indepdence sampler in that new proposed parameters are not dependent upon the current value of the chain.</p>
<p>After completion, the function returns convergence diagnostics (effective sample size, acceptance proportion ratio), MCMC chain, and the median value and 95% Highest Posterior Density of the estimated parameter.</p>
<p>The function ‘transformPhylo.MCMC’ allows for the estimation of model parameters using Bayesian statistics. Models of lambda, delta, kappa, OU, ACDC, psi, and multi-psi can currently be modelled using transformPhylo.MCMC. Additionally, Pagel’s lambda can be optimised alongside parameters and nested modes in the same way as <code>transformPhylo.ML</code>.</p>
<p>We will run an MCMC chain of 1000 generations to estimate Pagel’s lambda and discarding the first 10% (‘200 generations (’burn.in = 0.1’). All the models use a ‘uniform’ prior for each of the parameters. For lambda, this is a uniform distribution between 0 and 1 (although lambda can reach slightly higher than one), meaning we think all potential values are equally likely. To obtain identical results wel will set ‘random.start=FALSE’, if this is set to TRUE a random start value is taken from the system time</p>
<pre class="r"><code>set.seed(12)  # set seed so run will be identical - for example use only
lambda.mcmc &lt;- transformPhylo.MCMC(y = male.length.clade, phy = phy.clade, 
    model = &quot;lambda&quot;, mcmc.iteration = 2000, burn.in = 0.25, 
    random.start = FALSE, sample.every = 1)</code></pre>
<p>We can know check the posterior estimate of lambda and convergence of the model. The median and 95 Highest Posterior Density (HPD) is output by the model. Some diagnostics are output as standard: Effective Sample Size (ESS) and acceptance rate. We aim for an ESS of at least 200 and an acceptance rate around 0.44</p>
<pre class="r"><code>lambda.mcmc[1:4]</code></pre>
<pre><code>## $median
##    Lambda 
## 0.7893048 
## 
## $`95.HPD`
## lower 95% HPD upper 95% HPD 
##     0.5169347     0.9649743 
## 
## $ESS
##   Lambda 
## 356.8689 
## 
## $acceptance.rate
## [1] 0.3544304</code></pre>
<p>Our lambda median value is 0.79 but there is a large 95% HPD (0.52-0.96). The ESS and acceptance rate look ok. We can also plot the trace from the MCMC chain - this could look better - running for more generations would help</p>
<pre class="r"><code>mcmc.plot(lambda.mcmc)</code></pre>
<div class="figure">
<img src="/vignettes/figures/plot12-1.png" alt="Figure 12. MCMC trace for Pagel&#39;s lambda" width="1000" />
<p class="caption">
Figure 12. MCMC trace for Pagel’s lambda
</p>
</div>
</div>
<div id="character-displacement-models" class="section level1">
<h1>Character displacement models</h1>
<p>Magnus Clarke <em>et al.</em> (2017) introduced a character displacement model in which inter-specific competition can drive trait change. This model estimates a parameter ‘a’ that drives the strength of inter-specific competition, alongside a Brownian motion model with parameter estimation of the trait variance. If a=0 the model is equivalent to Brownian motion, and larger values of a drive trait evolution away from the values of inter-specific competitors.</p>
<p>The character displacement model employs an approximate Bayesian computation (ABC) approach, in which many datasets are simulated based on the known tree using a range of parameter values for <em>a</em> and the trait variance. These simulations then are compared to the empirical data to estimate the ‘best-fitting’ parameters of the Brownian motion process variance, and the character displacement parameter <em>a</em>.</p>
<p>First data are simulated on the known tree, allowing for a range of variance (sigma) and <em>a</em> values with both sample from a uniform distribution between 0 and 8. For brevity, we will use 100 simulations only. For actual analyses, many more iterations would be required, perhaps 1 million (Clarke <em>et al</em> 2017). Note this process can be made parallel on Mac and Linux systems by using the ‘mc.cores’ argument, but here we will use one core only.</p>
<pre class="r"><code>data(finches)
emp.tree &lt;- finch.tree
emp.data &lt;- finch.data
param.simulation &lt;- chr.disp.param(emp.tree, n.sim = 100, max.sigma = 8, 
    max.a = 8, ntraits = 1, mc.cores = 1)</code></pre>
<p>We can then compare these simulated data with the empirical data using the function ‘chr.disp.lrt’. We will use only 75 simulations from the posterior, this value can be guided by simulations (see Clarke et al. 2017)</p>
<pre class="r"><code>chr.disp.lrt(emp.tree = emp.tree, emp.data = emp.data, param.out = param.simulation, 
    posteriorSize = 75)</code></pre>
<pre><code>## $estimates
##       Brownian motion Character displacement model
## sigma        1.226667                     4.800000
## a            0.000000                     4.906667
## 
## $likelihood
##   log.h.0.lik. log.h.1.lik. likelihood.ratio.test   p.value
## 1    -5.179175    -3.706954              2.944443 0.9138266</code></pre>
<p>The output shows the ‘estimates’ for hypothesis 0 (Brownian motion) and hypothesis 1 (character displacement) with the variance and a values summarised (a is 0 in the Brownian motion model, by definition). The second list element ‘likelihood.ratio.test’ shows the likelihood of each model, the value of the likelihood-ratio test statistic, and the <em>p</em> value (here the character displacement is not supported over the character displacement model).</p>
</div>
<div id="fast-estimation-of-phylogenetic-generalised-least-squares" class="section level1">
<h1>Fast estimation of Phylogenetic Generalised Least Squares</h1>
<p>The package <em>caper</em> (Orme <em>et al</em> 2018) offers an excellent model to run Phylogenetic Generalised Least Squares (PGLS) models, but these are based-upon Generalised Least Squares (using variance-covariance matrices) which are substantially slower than using indpendent contrasts (Freckleton 2012).</p>
<p>In motmot, code allows for continuous PGLS models can be estimated using contrasts - this gives identical results to <em>caper</em> but is substantially faster, as is shown below. At current only continuous data is allowed in the models for motmot, so if any of the input data are not continuous CAPER or similar should be used. Additionally motmot only estimates Pagel’s lambda rather than other models, such as Kappa as offered by CAPER</p>
<pre class="r"><code># Data and phylogeny
data(anolis.tree)
anolis.tree$node.label &lt;- NULL
set.seed(3492)
lm.data &lt;- transformPhylo.sim(phy = anolis.tree, n = 2, model = &quot;bm&quot;)
dat &lt;- data.frame(x = lm.data[, 1], y = lm.data[, 2], names = anolis.tree$tip, 
    row.names = anolis.tree$tip)
# pgls from CAPER with matrix inversion
library(caper)</code></pre>
<pre><code>## Loading required package: MASS</code></pre>
<pre><code>## Loading required package: mvtnorm</code></pre>
<pre class="r"><code>comp.dat &lt;- comparative.data(anolis.tree, dat, names)
time.now &lt;- Sys.time()
matrix.inv.caper &lt;- pgls(y ~ x, data = comp.dat, lambda = &quot;ML&quot;)
pgls.time &lt;- Sys.time() - time.now
pgls.time</code></pre>
<pre><code>## Time difference of 0.9490931 secs</code></pre>
<pre class="r"><code>time.now &lt;- Sys.time()
picModel &lt;- pic.pgls(formula = y ~ x, phy = anolis.tree, y = dat, 
    lambda = &quot;ML&quot;, return.intercept.stat = FALSE)
pic.time &lt;- Sys.time() - time.now
pic.time</code></pre>
<pre><code>## Time difference of 0.2335591 secs</code></pre>
<p>The results are identical between the two methods</p>
<pre class="r"><code># from caper
summary(matrix.inv.caper)</code></pre>
<pre><code>## 
## Call:
## pgls(formula = y ~ x, data = comp.dat, lambda = &quot;ML&quot;)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -2.27281 -0.52137  0.08971  0.75199  2.21433 
## 
## Branch length transformations:
## 
## kappa  [Fix]  : 1.000
## lambda [ ML]  : 1.000
##    lower bound : 0.000, p = &lt; 2.22e-16
##    upper bound : 1.000, p = 1    
##    95.0% CI   : (0.873, NA)
## delta  [Fix]  : 1.000
## 
## Coefficients:
##              Estimate Std. Error t value Pr(&gt;|t|)
## (Intercept) -1.147535   2.685655 -0.4273   0.6697
## x            0.106595   0.073089  1.4584   0.1466
## 
## Residual standard error: 0.9448 on 163 degrees of freedom
## Multiple R-squared: 0.01288, Adjusted R-squared: 0.006825 
## F-statistic: 2.127 on 1 and 163 DF,  p-value: 0.1466</code></pre>
<pre class="r"><code># from MOTMOT
picModel</code></pre>
<pre><code>## $model
## 
## Call:
## lm(formula = formula, data = pic.data, x = TRUE, y = TRUE)
## 
## Coefficients:
##      x  
## 0.1066  
## 
## 
## $model.summary
## 
## Call:
## lm(formula = formula, data = pic.data, x = TRUE, y = TRUE)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.1581 -0.8276 -0.1225  0.6007  2.3315 
## 
## Coefficients:
##   Estimate Std. Error t value Pr(&gt;|t|)
## x  0.10660    0.07309   1.458    0.147
## 
## Residual standard error: 0.9448 on 163 degrees of freedom
## Multiple R-squared:  0.01288,    Adjusted R-squared:  0.006825 
## F-statistic: 2.127 on 1 and 163 DF,  p-value: 0.1466
## 
## 
## $intercept
## [1] -1.147535
## 
## $lambda
## [1] 1
## 
## $logLikelihood
##           [,1]
## [1,] -523.2653
## 
## $AIC
##          [,1]
## [1,] 1047.531</code></pre>
<p><strong>References</strong></p>
<ul>
<li>Blomberg SP, Garland T, and Ives AR. 2003. Testing for phylogenetic signal in comparative data: behavorial traits more labile. <em>Evolution</em> 57, 717–45.</li>
<li>Butler MA, and King AA. 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. <em>The American Naturalist</em> 164, 683-695.</li>
<li>Cavalli‐Sforza, LL, and Edwards AWF. 1967 Phylogenetic analysis: models and estimation procedures. <em>Evolution</em> 21, 550-570.</li>
<li>Clarke M, Thomas GH, and Freckleton RP. 2017. Trait evolution in adaptive radiations: modeling and measuring interspecific competition on phylogenies. <em>The American Naturalist</em> 189, 121-137.</li>
<li>Cooper N, Thomas GH, Venditti C, Meade A, &amp; Freckleton RP. 2016. A cautionary note on the use of Ornstein Uhlenbeck models in macroevolutionary studies. <em>Biological Journal of the Linnean Society</em> 118, 64-77.</li>
<li>Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters.]<em>American journal of human genetics</em> 25, 471.</li>
<li>Felsenstein J. 1985. Phylogenies and the comparative method.<em>The American Naturalist</em> 125, 1-15.</li>
<li>Freckleton RP. 2012. Fast likelihood calculations for comparative analyses. <em>Methods in Ecology and Evolution</em> 3, 940-947.</li>
<li>Hansen TF, 1997. Stabilizing selection and the comparative analysis of adaptation. <em>Evolution</em> 51, 1341-1351.</li>
<li>Harmon LJ, <em>et al.</em> 2010. Early bursts of body size and shape evolution are rare in comparative data. <em>Evolution</em> 64, 2385–96.</li>
<li>Ingram T. 2011. Speciation along a depth gradient in a marine adaptive radiation. <em>Proceedings of the Royal Society of London B: Biological Sciences</em> 278, 613-618.</li>
<li>Ingram T <em>et al</em>. 2016. Comparative tests of the role of dewlap size in <em>Anolis</em> lizard speciation. <em>Proceedings of the Royal Society of London B: Biological Sciences</em>, 283, 20162199.</li>
<li>O’Meara BC, Ané C, Sanderson MJ, and Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. <em>Evolution</em> 60, 922–33.</li>
<li>Orme D, Freckleton RP, Thomas GH, Petzoldt T, Fritz S, Isaac N, and Pearse W. 2018. caper: Comparative Analyses of Phylogenetics and Evolution in R. R package version 1.0.1.</li>
<li>Paradis E, Schliep K, and Schwartz R. 2018. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. <em>Bioinformatics</em>.</li>
<li>Pagel, M. Inferring evolutionary processes from phylogenies. 1997. <em>Zoologica Scripta</em> 26, 331-348.</li>
<li>Pagel, M. 1999. Inferring the historical patterns of biological evolution. <em>Nature</em> 401, 877.</li>
<li>Puttick, MN. 2018. Mixed evidence for early bursts of morphological evolution in extant clades. <em>Journal of Evolutionary Biology</em> 31, 502-515.</li>
<li>Slater GJ. 2013. Phylogenetic evidence for a shift in the mode of mammalian body size evolution at the Cretaceous‐Palaeogene boundary. Methods in Ecology and Evolution, 4, 734-744.</li>
<li>Thomas GH, and Freckleton RP. 2012. MOTMOT: Models of trait macroevolution on trees. <em>Methods in Ecology and Evolution</em> 3, 145–51.</li>
<li>Thomas GH, Freckleton RP, and Székely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. <em>Proceedings of the Royal Society B: Biological Sciences</em> 273, 1619–24.</li>
</ul>
</div>
</div>
