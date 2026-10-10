# Literature review for `cnavier_schemes.tex`

A structured literature search, not a systematic review in the PRISMA sense. Search date: 2026-10-09. The search was done in three strands, run in parallel, each logged below. Every candidate was checked against Crossref, arXiv or publisher metadata, and its abstract (or full text, where marked) was read before inclusion. Claims in the paper about a cited work are limited to what was verified at the level stated in the table.

**Strands**

- **A.** Comparisons of compact or high-order explicit finite differences with pseudospectral methods in turbulence. Forward citations of San & Staples 2012 and Laizet & Lamballais 2009, plus keyword searches. About 1,050 citation records and 600 keyword hits screened by title, about 60 abstracts read.
- **B.** Aliasing versus truncation, and energy pile-up or thermalization at a spectral cutoff. Forward citations of Kravchenko & Moin 1997, Frisch et al. 2008, Cichowlas et al. 2005 and Fox & Orszag 1973, plus keyword searches. About 990 records screened, about 60 abstracts read.
- **C.** Resolution criteria and accuracy-per-cost comparisons of schemes. About 600 records screened, 60 examined in detail.

Sources: OpenAlex (forward citations, search, abstracts), Semantic Scholar (citations, abstracts; rate-limited for strand C), Crossref (field checks), the arXiv API and PDFs, WebSearch. ScienceDirect, SpringerLink and Google Scholar mostly refused automated access, which is why some items are verified only to the abstract or metadata.

## What the review changed in the paper

1. **Resolution axis.** The classical two-dimensional criterion, k_max·l_η with l_η = (ν³/η)^(1/6) (Lunasin et al. 2007, full text checked; η = ν⟨|∇ω|²⟩ as in our code; their table gives k_max = √2N/3 in units of 2π on a unit box, and we infer that the product uses angular wavenumbers, as the 3-D k_max·η does, since shell units would put ~38 grid cells per l_η in their resolved runs), orders the ranking of schemes monotonically in our ten cases, with boundaries near its thresholds of 1 (resolved) and 2 (well resolved). It replaces the cutoff energy χ as the paper's main axis. χ orders the cases only approximately. A cutoff-to-peak spectral ratio is also already used as a pass/fail convergence diagnostic (Baty 2026, following García Morillo & Alexakis 2025).
2. **Accuracy metric.** The resolved bandwidth has a precedent in the "effective spectral bandwidth" of Kritsuk et al. 2011 (full text of §5.4 checked: ±25 % of a filtered reference, with cost left out). The paper cites it and claims only the tighter tolerances, the cross-validated references and the combination with cost.
3. **Cost metric.** Kawai et al. 2026 combine effective resolution from spectra with resource use, including a stability-limited time step (abstract only; the full text was behind a login). The paper cites it as precedent. Raeth & Hallatschek 2024 show the time-step limit depends on the form of the nonlinear term, so the paper states that its ratios hold for the skew-symmetric form.
4. **Pile-up.** No paper found runs the same marginally resolved flow with and without 3/2 padding, so the control appears to be new. The truncation interpretation is not new: Ishihara et al. 2018 (full-text passage), Bardos & Tadmor 2015, and the Galerkin-truncation literature (Cichowlas et al. 2005; Frisch et al. 2008; Ray et al. 2011). The conventional aliasing view is stated in Fontana et al. 2020. The paper restricts its conclusion to the skew-symmetric form (Zang 1991; Blaisdell et al. 1996).
5. **Two-dimensional comparisons.** Fox & Orszag 1973, Herring et al. 1974 and Browning & Kreiss 1989 were added as the earlier two-dimensional spectral versus finite-difference comparisons. Vreman et al. 1996 was added for the reversal of scheme ranking with resolution.

## Cited, with verification level

| Key | Work | Verified | Used for |
|---|---|---|---|
| san2012 | San & Staples, Comput. Fluids 63 (2012) 105 | full text (arXiv:1212.0920) | the work the paper responds to |
| lunasin2007 | Lunasin, Kurien, Taylor & Titi, J. Turbul. 8 (2007) N30 | full-text passage, Eqs. (29)–(30) | k_max·l_η and its thresholds |
| kritsuk2011 | Kritsuk et al., ApJ 737 (2011) 13 | full text §5.4, §6 | precedent for the bandwidth metric |
| ishihara2018 | Ishihara et al., ApJ 854 (2018) 81 | full-text passage §2 | truncation reading of the pile-up |
| fontana2020 | Fontana, Bruno, Mininni & Dmitruk, CPC 256 (2020) 107482 | full-text passage | the aliasing view |
| baty2026 | Baty, arXiv:2604.02065 (2026), following García Morillo & Alexakis, JFM 1007 (2025) R3 | Baty full text §3; García Morillo & Alexakis abstract | cutoff-to-peak ratio as a diagnostic |
| lambert2026 | Lambert, Reneuve & Augier, arXiv:2603.08892 (2026) | abstract, first half of the text | phase-shift dealiasing cost |
| bardos2015 | Bardos & Tadmor, Numer. Math. 129 (2015) 749 | abstract | 2/3 method shares the high-mode pathology |
| cichowlas2005, ray2011, frisch2008 | PRL 95 (2005) 264502; PRE 84 (2011) 016301; PRL 101 (2008) 144501 | abstract | Galerkin-truncation thermalization |
| herring1974 | Herring, Orszag, Kraichnan & Fox, JFM 66 (1974) 417 | abstract | 2-D spectral vs grid-point accuracy |
| browning1989 | Browning & Kreiss, Math. Comp. 52 (1989) 369 | abstract | 2-D converged references, scheme comparison |
| vreman1996 | Vreman, Geurts & Kuerten, IJNMF 22 (1996) 297 | abstract | ranking reverses with resolution |
| yeung2018 | Yeung, Sreenivasan & Pope, PRFluids 3 (2018) 064603 | abstract | k_max·η in 3-D |
| kawai2026 | Kawai et al., J. Meteorol. Soc. Japan 104 (2026) 20 | abstract only | precedent for the cost metric |
| raeth2024 | Raeth & Hallatschek, Phys. Fluids 36 (2024) 105167 | abstract | time step depends on the nonlinear form |
| sinhababu2021 | Sinhababu & Ayyalasomayajula, Math. Comput. Simul. 182 (2021) 116 | abstract | 2/3 costs more than 3/2 at equal accuracy |
| hou2007 | Hou & Li, JCP 226 (2007) 379 | abstract | high-order Fourier filtering |
| fornberg1987 | Fornberg, Geophysics 52 (1987) 483 | abstract | early cost comparison |
| fox1973 | Fox & Orszag, JCP 11 (1973) 612 | metadata (abstract withheld) | early 2-D pseudospectral turbulence |
| zang1991, blaisdell1996 | Appl. Numer. Math. 7 (1991) 27; 21 (1996) 207 | metadata (abstracts withheld) | skew-symmetric form reduces aliasing; **read before submission** |
| fedioun2001, park2004, yalla2021, yamamoto2026, akkurt2025, canuto2006, rodhiya2026 | see the paper | checked in the earlier pass (PR #55) | |

## Considered and not cited

- **Lee & Seo 2002**, JCP 183, 438: no abstract or text was accessible. Citing papers describe it as a windowed sinc (Fourier) derivative stencil.
- **Gotoh, Hatanaka & Miura 2012**, JCP 231, 7398 (hybrid Fourier/combined compact code, 26–77 % cheaper at spectral accuracy): only a search snippet was available.
- **Sekimoto, Dong & Jiménez 2016**: effective cutoff of a compact scheme from its modified wavenumber. Peripheral.
- **Pope 2000** (k_max·η ≥ 1.5): confirmed only through citing papers; Yeung et al. 2018 is cited instead.
- **Bustamante & Brachet 2012** and **Sulem, Sulem & Frisch 1983** (analyticity-strip criterion): related resolution checks read from the spectral tail, but inviscid and single-method.
- **Pirozzoli 2007**, **Kreiss & Oliger 1972**, **Johnsen et al. 2010**, **Skamarock 2004**, **Kent et al. 2014**, **Capuano et al. 2023**: cost or effective-resolution methodology in other settings (wave propagation, NWP, compressible shocks, sphere flow). Some were reachable only through secondary sources.
- **Mostipan & Troshin 2026** and **Song et al. 2024**: explicit finite-difference dispersion and a compressible compact framework. Peripheral.
- **Margairaz et al. 2018**, **Bowman & Roberts 2011**, **Bos & Bertoglio 2006**, **Agrawal et al. 2020**, **Yamazaki, Ishihara & Kaneda 2002**: dealiasing implementations and bottleneck closures. Related but not needed for the paper's claims.

Note for double-blind review: the cnavier repository's issues and pull requests are indexed by web search for these topics.

## Search logs

### Strand A: compact / high-order FD vs pseudospectral

| # | Source / call | Hits returned | Screened |
|---|---|---|---|
| 1 | OpenAlex `works/doi:` for San & Staples 2012, Laizet & Lamballais 2009, Lee & Seo 2002. Resolved to W1994203821, W2120318017, W2047436969. All three have `abstract_inverted_index = null`. | 3 | 3 |
| 2 | OpenAlex `filter=cites:W1994203821` (San & Staples forward citations) | 63 | 63 (all titles; abstracts for candidates) |
| 3 | Semantic Scholar `paper/DOI:10.1016/j.compfluid.2012.04.006/citations` | 70 | 70 (12 not in OpenAlex, screened separately) |
| 4 | OpenAlex `filter=cites:W2120318017` (Laizet & Lamballais forward citations), pages 1–3 | 410 | 410 titles (keyword filter followed by a manual pass of the full list) |
| 5 | Semantic Scholar `paper/DOI:10.1016/j.jcp.2009.05.010/citations` | 443 | ~105 not in OpenAlex screened; the rest duplicated #4 |
| 6 | OpenAlex `filter=cites:W2047436969` (Lee & Seo forward citations) | 44 | 44 |
| 7 | Semantic Scholar Lee & Seo record, plus citations with `contexts` | 1 record; ~25 citing contexts | all |
| 8 | Crossref `works/10.1006/jcph.2002.7201`; OpenAlex referenced_works of Lee & Seo (22 refs resolved) | 1 / 22 | all |
| 9 | WebFetch ScienceDirect PII S0021999102972010 and ResearchGate page for Lee & Seo | 403 both | – |
| 10 | Google Scholar (curl) and CORE API search for Lee & Seo | blocked / redirect error | – |
| 11 | WebSearch "A new compact spectral scheme for turbulence simulations" Lee Seo | 10 | 10 |
| 12 | WebSearch "Lee and Seo" 2002 compact spectral scheme … | ~28 | 28 |
| 13 | WebSearch Changhoon Lee Youngchwa Seo compact spectral scheme sinc window … | ~28 | 28 |
| 14 | OpenAlex search "compact finite difference pseudospectral comparison turbulence" | 629 | top 40 |
| 15 | OpenAlex search "compact scheme spectral method under-resolved turbulence" | 8903 | top 40 |
| 16 | OpenAlex search "two-dimensional turbulence high-order finite difference spectral comparison" | 33293 | top 40 |
| 17 | OpenAlex search "spectral-like resolution compact scheme direct numerical simulation accuracy cost" | 15351 | top 40 |
| 18 | OpenAlex search "Arakawa compact pseudospectral two-dimensional turbulence" | 23 | 23 |
| 19 | OpenAlex search "finite difference versus spectral isotropic turbulence resolution" | 10409 | top 40 |
| 20 | OpenAlex search "compact difference Fourier spectral homogeneous turbulence comparison accuracy" | 2794 | top 40 |
| 21 | OpenAlex search "aliasing error finite difference spectral turbulence simulation" | 2692 | top 40 |
| 22 | OpenAlex search "compact spectral scheme turbulence Lee Seo" | – | top 10 |
| 23 | arXiv API: `abs:compact AND abs:pseudospectral AND abs:turbulence` | 2 | 2 |
| 24 | arXiv API: `abs:"compact scheme" AND abs:spectral AND abs:"under-resolved"` | 1 | 1 |
| 25 | arXiv API: `abs:"finite difference" AND abs:pseudo-spectral AND abs:turbulence AND abs:accuracy` | 1 | 1 |
| 26 | arXiv API: `abs:"two-dimensional turbulence" AND abs:"finite difference" AND abs:spectral` | 0 | – |
| 27 | arXiv API: `abs:"modified wavenumber" AND abs:turbulence` | 0 | – |
| 28 | arXiv API: `abs:dealiasing AND abs:"finite difference" AND abs:turbulence` | 2 | 2 |
| 29 | arXiv API: `abs:compact AND abs:"finite difference" AND abs:spectral AND abs:turbulence` | 8 | 8 |
| 30 | arXiv API: `abs:"spectral-like" AND abs:turbulence` | 7 | 7 |
| 31 | arXiv API: `abs:"finite-difference" AND abs:"pseudo-spectral" AND abs:turbulence` | 6 | 6 |
| 32 | arXiv API: `abs:"finite difference" AND abs:"pseudospectral" AND abs:turbulence` | 2 | 2 |
| 33 | arXiv API: `abs:"cell Reynolds number" AND abs:turbulence` | 1 (San & Staples only) | 1 |
| 34 | arXiv API: `abs:"effective resolution" AND abs:turbulence AND abs:spectral` | 1 | 1 |
| 35 | arXiv API: `abs:"compact spectral scheme"`; `ti:"shortcomings" AND au:Tadmor` | 0 / 1 | 1 |
| 36 | WebSearch compact FD vs pseudo-spectral DNS isotropic turbulence comparison energy spectrum | 9 | 9 |
| 37 | WebSearch Gotoh Hatanaka Miura spectral compact difference hybrid | 9 | 9 |
| 38 | WebSearch 2-D decaying turbulence comparison FD Arakawa pseudospectral | 10 | 10 |
| 39 | WebSearch "Accuracy and computational efficiency of dealiasing schemes …" | 9 | 9 |
| 40 | WebSearch Herring Orszag Kraichnan Fox 1974 … | 10 | 10 |
| 41 | WebSearch cost vs accuracy high-order FD / compact / pseudospectral GPU | ~18 | 18 |
| 42 | WebSearch spectral blocking / pile-up near cutoff in dealiased pseudospectral 2-D | ~27 | 27 |
| 43 | WebSearch vorticity–streamfunction 2-D sixth-order compact vs Fourier | 9 | 9 |
| 44 | WebSearch resolution criterion k_max·η for high-order FD vs spectral | ~28 | 28 |
| 45 | WebSearch GPU 2-D turbulence pseudospectral vs FD wall-clock | ~18 | 18 |
| 46 | WebSearch Kreiss & Oliger 1972 | ~29 | 29 |
| 47 | WebSearch resolution indicator: spectrum at cutoff / peak ratio | ~27 | 27 |
| 48 | Crossref `query.bibliographic` DOI lookups (12 queries) and `works/DOI` checks (~50 DOIs); Semantic Scholar batch abstracts (16 DOIs) | – | – |
| 49 | Full-text checks (arXiv PDFs): Baty 2604.02065, Sekimoto et al. 1601.01646, Laizet et al. 1409.3621 | 3 | 3 |

Overall: about 1,050 forward-citation records and about 600 keyword hits were screened by title. Roughly 60 abstracts were read.

### Strand B: aliasing vs truncation, pile-up

| # | Source / call | Query or filter | Hits screened |
|---|---|---|---|
| 1 | OpenAlex `works/doi:10.1006/jcph.1996.5597` | Resolved the Kravchenko & Moin 1997 ID to W2034644845 | 1 |
| 2 | OpenAlex `works?filter=cites:W2034644845` (cursor paging, all pages) | All forward citations of Kravchenko & Moin 1997 | 534 titles. Keyword pre-filter (alias, truncat, pile, thermali, cutoff, pseudo-spectral, skew, 3/2, 2/3, under-resolved, bottleneck, compact, Fourier, two-dimensional) left 123, which were read by title. About 12 abstracts were read. |
| 3 | OpenAlex `works/doi:10.1103/PhysRevLett.101.144501` | Resolved the Frisch et al. 2008 ID to W2037896057 | 1 |
| 4 | OpenAlex `works?filter=cites:W2037896057` (all pages) | All forward citations of Frisch et al. 2008 | 190 titles, all read by title. About 15 abstracts were read. |
| 5 | OpenAlex `works?filter=cites:W2064882198` (all pages) | Forward citations of Cichowlas et al. 2005, screened for 2-D, aliasing, under-resolution and numerics | 170 titles. A pre-filter left about 100, all read by title. |
| 6 | OpenAlex `works?filter=cites:W2093328087` (all pages) | Forward citations of Fox & Orszag 1973 (2-D aliasing) | 93 titles, all read by title |
| 7 | OpenAlex `works?search=` | "Galerkin truncation thermalization two-dimensional turbulence" | 20 |
| 8 | OpenAlex search | "energy pile-up spectral cutoff under-resolved direct numerical simulation" | 20 |
| 9 | OpenAlex search | "3/2 rule 2/3 rule dealiasing comparison" | 20 |
| 10 | OpenAlex search | "aliasing error truncation error spectral method turbulence" | 20 |
| 11 | OpenAlex search | "skew-symmetric form aliasing error reduction turbulence" | 20 |
| 12 | OpenAlex search | "bottleneck effect spectral truncation thermalization" | 20 |
| 13 | OpenAlex search | "two-dimensional turbulence enstrophy pile-up truncation wavenumber" | 15 |
| 14 | OpenAlex search | "phase shift dealiasing pseudospectral cost" | 15 |
| 15 | OpenAlex search | "spectral filter versus dealiasing turbulence simulation" | 15 |
| 16 | OpenAlex search | "effects of finite spatial resolution direct numerical simulation turbulence spectrum high wavenumber" | 15 |
| 17 | OpenAlex search | "spectrally truncated inviscid turbulence" | 15 |
| 18 | OpenAlex search | "dealiasing pseudospectral two-dimensional turbulence" | 15 |
| 19 | OpenAlex search | "spectral blocking"; "spectral blocking dealiasing aliasing"; "effects of finite spatial and temporal resolution …"; "under-resolved spectral simulation energy accumulation highest wavenumbers dealiased" | 60 |
| 20 | WebSearch | `"spectral blocking" Boyd aliasing 2/3 rule does not cure energy pile-up highest wavenumbers` | about 30 links. Led to Boyd's book (full PDF read, ch. 11). |
| 21 | WebSearch | `Rodhiya Bhattacharya Verma "Relative accuracy …"` | 9 |
| 22 | WebSearch | `Yeung Sreenivasan Pope 2018 "finite spatial and temporal resolution" spectrum upturn near kmax …` | about 30. Led to Ishihara et al. 2018. |
| 23 | WebSearch | `"pile-up" "caused by the wavenumber truncation" Fourier spectral method "finite difference"` | about 20. Confirmed Ishihara et al. 2018 and found Fontana et al. 2020 and Alhawwary & Wang 2018. |
| 24 | WebSearch | `Fourier spectral DNS energy spectrum "pile-up" near kmax truncation dealiased Ishihara Kaneda …` | about 10. Led to Yamazaki et al. 2002 and Kaneda & Ishihara 2006. |
| 25 | WebSearch | `two-dimensional turbulence pseudospectral under-resolved enstrophy accumulation at truncation wavenumber …` | about 20 |
| 26 | WebSearch | `aliased versus dealiased pseudospectral simulation two-dimensional turbulence … "3/2 rule" …` | about 25. Led to Fox & Orszag 1973 and Bowman's talk slides. |
| 27 | arXiv API `id_list=` | 2603.08892, 1308.5314, 2209.05046, 2512.10788, 2508.10808, 2109.00255, 1008.1366, 2402.17688, 1801.08805, 1802.02719, 2004.06274, 0710.4100, 1409.3621, nlin/0502036, 2509.18512, 2301.09900 | 16 records |
| 28 | Crossref `works/{DOI}` and `query.bibliographic` | Field check for every included item. Also used to recover the DOIs of Bardos & Tadmor 2015, Kaneda & Ishihara 2006, Maltrud & Vallis 1993, Fox & Orszag 1973, Blaisdell et al. 1996, Zang 1991 and Patterson & Orszag 1971. | about 35 |
| 29 | Semantic Scholar `paper/DOI:` | Abstracts where OpenAlex had none: Hou & Li 2007, Bardos & Tadmor 2015, Yeung et al. 2018. The publisher withheld the abstract for Blaisdell 1996, Zang 1991, Yamazaki 2002, Fox & Orszag 1973 and Maltrud & Vallis 1993. | 9 |
| 30 | Full-text reads (PDF or HTML) | Boyd 2001 ch. 11; Ishihara et al. 2018 (arXiv 1801.08805); Fontana et al. 2020 (arXiv 2002.01392); Lambert et al. 2026 (arXiv 2603.08892, first half); Rodhiya et al. (arXiv 2508.10808, via summarising fetch); Bowman talk slides "How important is dealiasing"; Seshasayanan arXiv 2301.09900 (checked, not relevant) | 7 |

In total, about 990 records were screened by title and about 60 abstracts were read.

### Strand C: resolution criteria and cost

| # | Query / call | Source | Hits screened | Notes |
|---|---|---|---|---|
| 1 | `effective resolution kinetic energy spectrum numerical scheme` | OpenAlex `works?search=` | 15 of 38,406 | Found Soufflet 2016. Most hits were atmospheric or ocean spectra papers. |
| 2 | `resolution requirement direct numerical simulation k_max eta` | OpenAlex | 15 of 31 | Mostly noise. The first call was rate-limited and was retried. |
| 3 | `cell Reynolds number resolution criterion turbulence` | OpenAlex | 15 of 12,659 | Noise. Nothing relevant. |
| 4 | `resolving efficiency finite difference scheme cost` | OpenAlex | 15 of 84,342 | Noise. Nothing relevant. |
| 5 | `cost accuracy comparison compact scheme spectral GPU` | OpenAlex | 15 of 8,481 | Found Capuano et al. 2023 (cost vs accuracy). |
| 6 | `two-dimensional turbulence resolution enstrophy dissipation scale grid` | OpenAlex | 15 of 2,400 | Found Thuburn et al. 2014. |
| 7 | `effective resolution numerical scheme turbulence spectrum dissipation range` | OpenAlex | 15 of 9,250 | Found Kritsuk 2011 (ApJ), Maulik & San 2017 and astrophysics bottleneck papers. |
| 8 | `spatial resolution DNS small scale statistics extreme events` | OpenAlex | 15 of 2,064 | Background only. |
| 9 | title_and_abstract filter: `"effective resolution" turbulence scheme spectra` | OpenAlex | 7 of 7 | Found **Kawai et al. 2026 (JMSJ)**. |
| 10 | title_and_abstract filter: `"spectral bandwidth" turbulence numerical methods comparison` | OpenAlex | 2 of 2 | **Kritsuk et al. 2011** |
| 11 | title_and_abstract filter: `"resolution requirements" "direct numerical simulation" turbulence spectral` | OpenAlex | 12 of 13 | Peng et al. 2010 (LBM vs pseudospectral) |
| 12 | title_and_abstract filter: `"cell Reynolds number" turbulence high-order` | OpenAlex | 3 of 3 | Only San & Staples 2012 (already cited) and a 2009 compressible upwind paper |
| 13 | title_and_abstract filter: `"resolving efficiency"` | OpenAlex | 12 of 86 | Mostly noise. Optimised WENO and flux-reconstruction papers. |
| 14 | title_and_abstract filter: `compact pseudospectral efficiency turbulence "computational cost"` | OpenAlex | 0 | none |
| 15 | five keyword queries (cutoff energy criterion, effective resolution FD vs spectral, 2-D resolution, compact vs pseudospectral efficiency, cell Reynolds number) | Semantic Scholar `paper/search` | 0 | HTTP 429 every time (rate limit, no key). Replaced by the OpenAlex filters (rows 9 to 14) and WebSearch. |
| 16 | Yeung Sreenivasan Pope 2018 PRFluids finite resolution | WebSearch | 10 | Confirmed the DOI. Pointed to Donzis 2008. |
| 17 | Donzis Yeung Sreenivasan 2008 resolution effects | WebSearch | 10 | |
| 18 | Watanabe Gotoh 2007 inertial-range intermittency accuracy DNS | WebSearch | 10 | |
| 19 | Skamarock 2004 KE spectra effective resolution | WebSearch | 9 | Also found Klaver 2020 and the AROME/Meso-NH paper. |
| 20 | analyticity strip width spectral resolution criterion | WebSearch | 19 | Found Sulem, Sulem & Frisch 1983 and Bustamante & Brachet 2012. Also surfaced Rodhiya et al. (already cited). |
| 21 | Pirozzoli 2007 resolving efficiency cost | WebSearch | 28 | Abstract not reachable (ScienceDirect returned 403). Content is secondary. |
| 22 | Moura Sherwin Peiró 2015 "1% rule" | WebSearch | 10 | |
| 23 | Ghosal 1996 numerical errors LES spectrum cutoff | WebSearch | 10 | Also Chow & Moin 2003 and Park & Mahesh 2007. |
| 24 | 2-D turbulence resolution requirement k_max, enstrophy dissipation wavenumber | WebSearch | 18 | Found Lunasin et al. 2007 (k_max·l_η criterion) and Boffetta & Ecke 2012. |
| 25 | `"spectral occupancy" OR "energy at the cutoff"` resolution criterion | WebSearch | 19 | No paper uses either phrase as a named criterion. Found Yamamoto & Tsuji arXiv:2605.22107 (99% energy criterion). |
| 26 | Sytine et al. 2000 PPM convergence | WebSearch | 20 | |
| 27 | pseudospectral vs compact cost per resolved wavenumber, time-step limit | WebSearch | 19 | Found Baltzer & Livescu 2020. Also surfaced cnavier issue #51 (the author's own repository; excluded). |
| 28 | Soufflet 2016; Johnsen 2010; Vermeire 2017; Kritsuk 2011 | WebSearch ×4 | about 40 | |
| 29 | Johnsen et al. 2010 full text | UNL DigitalCommons PDF | 1 | Got an HTML/bot page, not the PDF. The paper stays at abstract/secondary level. |
| 30 | Lele, Larsson, Bhagatwala & Moin 2009 ASP 406 | WebSearch | 20 | Title identified, page number not confirmed |
| 31 | Capuano et al. 2023 | WebSearch, then NSF-PAR PDF | 9, 1 | Full text read (cost metric) |
| 32 | GPU pseudospectral vs finite-difference performance | WebSearch | 9 | Found Biswas & Ganesh 2024 and Ravikumar et al. SC19 |
| 33 | Bracco & McWilliams 2010 | WebSearch | 20 | |
| 34 | Kreiss & Oliger 1972 | WebSearch | 20 | Also Fornberg 1987 |
| 35 | Frehlich & Sharman 2008; Schranner et al. 2015; Thuburn et al. 2014 | WebSearch ×3 | about 30 | |
| 36 | `"energy at the Nyquist" OR "energy at the cutoff"` ratio-to-peak criterion for scheme choice | WebSearch | 10 | **The only hits stating this exact criterion are the cnavier PRs #52 and #54 (the author's own work).** Also found Fernandez et al. 2019. |
| 37 | ratio of spectrum at k_max to spectral peak, "decades", 2-D | WebSearch | 29 | No formal E(k_max)/E_peak rule found. Found Luo, Fang & Fang 2025. |
| 38 | effective resolution of compact FD vs spectral reference, GPU cost, 2024–2025 | WebSearch | 47 | Found Kent et al. 2014, Maron 2004 and Cao et al. 2019. No head-to-head GPU compact vs pseudospectral cost study at equal accuracy. |
| 39 | Kent et al. 2014 Parts I/II | WebSearch | 19 | Abstracts |
| 40 | Pope 2000 k_max η ≥ 1.5 | WebSearch | 19 | Chapter DOI found. Criterion confirmed only through citing papers. |
| A | arXiv API `id_list=` 2508.10808, 1802.02719, 2512.20676, 2004.06274, 2006.13492, 2605.22107, 1909.02907, 1409.3621, 1103.5525, 1212.0920, 1112.1571, 2402.05478, 2206.04729, 1804.09712, astro-ph/0402606, 1905.03007, physics/0702196 | arXiv | 17 | Metadata and abstracts |
| B | Crossref `works/{DOI}` for 28 DOIs, plus `query.bibliographic` for 9 more | Crossref | 37 | All the DOIs below resolve and their fields match |
| C | OpenAlex `works/doi:` `abstract_inverted_index` for 18 DOIs | OpenAlex | 18 | 8 abstracts recovered |
| D | Kritsuk et al. 2011 arXiv PDF, read with `pdftotext` and grep | arXiv | 1 | Full text: §5.4 "Effective Spectral Bandwidth" and the Discussion |
| E | Lunasin et al. 2007 arXiv PDF | arXiv | 1 | Full text: the 2-D resolution criterion at lines 636–655 |

Total screened: about 600 records, 60 examined in detail, 18 included.
