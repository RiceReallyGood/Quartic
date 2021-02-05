# Exact Diagonalization

## Hamiltonian

$$
H = -t\sum_{i, \sigma}\left(c_{i,\sigma}^{\dagger}c_{i-1,\sigma} + c_{i,\sigma}^{\dagger} c_{i+1, \sigma}\right) + U \sum_{i} n_{i,\uparrow} n_{i, \downarrow}
$$

Fourier transform of operators
$$
\left\{
\begin{aligned}
c_{k,\sigma} & =\frac{1}{\sqrt{N_{s}}}\sum_{i}c_{i,\sigma}e^{-I\frac{2\pi k}{N_{s}}i}\\
c_{k,\sigma}^{\dagger} & =\frac{1}{\sqrt{N_{s}}}\sum_{i}c_{i,\sigma}^{\dagger}e^{I\frac{2\pi k}{N_{s}}i}
\end{aligned}
\right.
\qquad
\left\{
\begin{aligned}c_{i,\sigma} & =\frac{1}{\sqrt{N_{s}}}\sum_{k}c_{k,\sigma}e^{I\frac{2\pi k}{N_{s}}i}\\
c_{i,\sigma}^{\dagger} & =\frac{1}{\sqrt{N_{s}}}\sum_{k}c_{k,\sigma}^{\dagger}e^{-I\frac{2\pi k}{N_{s}}i}
\end{aligned}
\right.
$$
First term:
$$
\begin{aligned} & -t\sum_{i,\sigma}\left(c_{i,\sigma}^{\dagger}c_{i-1,\sigma}+c_{i,\sigma}^{\dagger}c_{i+1,\sigma}\right)\\
= & -t\sum_{i,\sigma}\left(\frac{1}{N_{s}}\sum_{k,k'}c_{k,\sigma}^{\dagger}c_{k',\sigma}e^{-I\frac{2\pi k}{N_{s}}i}e^{I\frac{2\pi k'}{N_{s}}\left(i-1\right)}+\frac{1}{N_{s}}\sum_{k,k'}c_{k,\sigma}^{\dagger}c_{k',\sigma}e^{-I\frac{2\pi k}{N_{s}}i}e^{I\frac{2\pi k'}{N_{s}}\left(i+1\right)}\right)\\
= & -\frac{t}{N_{s}}\sum_{k,k',\sigma}c_{k,\sigma}^{\dagger}c_{k',\sigma}N_{s}\delta_{kk'}2\cos\left(\frac{2\pi k}{N_{s}}\right)\\
= & \sum_{k,\sigma}\epsilon_{k}c_{k,\sigma}^{\dagger}c_{k,\sigma}
\end{aligned}
$$
where $\epsilon_{k}=-2t\cos\left(\frac{2\pi k}{N_{s}}\right)$

Second term:
$$
\begin{aligned} & U\sum_{i}n_{i,\uparrow}n_{i,\downarrow}\\
= & U\sum_{i}\frac{1}{N_{s}^{2}}\sum_{k_{1},k_{2},k_{3},k_{4}}c_{k_{1},\uparrow}^{\dagger}c_{k_{2},\uparrow}c_{k_{3},\downarrow}^{\dagger}c_{k_{4},\downarrow}e^{-I\frac{2\pi\left(k_{1}-k_{2}+k_{3}-k_{4}\right)}{N_{s}}i}\\
= & \frac{U}{N_{s}}\sum_{k_{1},k_{2},k_{3},k_{4}}c_{k_{1},\uparrow}^{\dagger}c_{k_{2},\uparrow}c_{k_{3},\downarrow}^{\dagger}c_{k_{4},\downarrow}\delta_{k_{1}+k_{3},k_{2}+k_{4}}^{\left(N_{s}\right)}
\end{aligned}
$$
So
$$
H = \sum_{k,\sigma}\epsilon_{k}c_{k,\sigma}^{\dagger}c_{k,\sigma} + \frac{U}{N_{s}}\sum_{k_{1},k_{2},k_{3},k_{4}}c_{k_{1},\uparrow}^{\dagger}c_{k_{2},\uparrow}c_{k_{3},\downarrow}^{\dagger}c_{k_{4},\downarrow}\delta_{k_{1}+k_{3},k_{2}+k_{4}}^{\left(N_{s}\right)}
$$
