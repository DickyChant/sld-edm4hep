# SLD mini-DST banks → EDM4hep

Mapping version 2026-09-26. Target layout: the collections used by the ALEPH/OPAL legacy
conversions (podio 1.7, EDM4hep 1.0), so that the same reader works for all three
experiments. Lengths cm → mm, energies GeV, magnetic field 0.6 T.

Every convention below was either taken from the jazelle_reader source (bank field
order) or established empirically on the released sample and then verified by the
read-back checks in `scripts/validate.py`. Items that could not be fixed are listed in
section 4 and are preserved raw so that nothing is lost.

## 1. Collections

| SLD bank | EDM4hep | fields |
|---|---|---|
| IEVENTH | `EventHeader` | runNumber, eventNumber, timeStamp = `evttime` (Unix seconds), weight; `trigger`, `evttype`, `header_word` as event parameters |
| PHBM | event parameters `ecm, decm, pol, dpol`; `PrimaryVertex` | position = `pos`·10, covMatrix = `dpos`·100 (6 packed lower-triangle terms), type = primary. `ecm` is −999 where unmeasured, 0 in simulation. `pol` is the signed longitudinal electron-beam polarization |
| PHPSUM + PHPOINT | `ReconstructedParticles` | momentum, energy = `ptot` (massless convention), mass 0, charge, referencePoint = (x,y,z)·10; relations `tracks` (PHPOINT.phchrg_id) and `clusters` (PHPOINT.phklus_id). PHPOINT has one row per PHPSUM row; a pointer of 0 means none |
| PHCHRG | `ChargedTracks` | type = `status` word (raw), chi2 = `chi2dt`, ndf = `ndfdt`, Nholes = `nmisht`, subdetectorHitNumbers = [`nhit`, `nhitv`, `nhite`, `nhitp`, `nmisht`, `nwrght`]; one `TrackState` (AtIP), see §2.1 |
| PHCHRG (extra) | `ChargedTracksDedxRaw` (int32), `ChargedTracksImcRaw` (int32), `ChargedTracksImpact2D/3D` (mm), `ChargedTracksImpact2DSignificance/3DSignificance` | raw `dedx` word, raw `imc`, `bnorm`, `b3norm`, `impact`, `impact3` |
| PHKLUS | `Clusters` | type = `status` (raw), energy = `eraw`, position = point on the nominal LAC front face along (cosθ, φ) (§2.2), iTheta, iPhi, shapeParameters = [cth, wcth, phi, wphi, nhit2, cth2, wcth2, phi2, whphi2, nhit3, cth3, wcth3, phi3, wphi3], subdetectorEnergies = `elayer[8]` |
| PHKELID | `ElectronID` (ParticleID, algorithmType 1, PDG 11) | type = `idstat`, likelihood = `prob`, parameters = [qp, em1, em12, had1, isolat, phiwid, thewid, em1x1, em2x2a, em2x2b, em3x3a, em3x3b]; `particle` → RP via PHPOINT.phkelid_id |
| PHWIC | `MuonID` (ParticleID, algorithmType 2, PDG 13) | type = `idstat`, likelihood = `missprob`, parameters = [nhit, nhit45, npat, nhitpat, nlayexp, nlaybey, chi2, ndf, matchChi2, matchNdf, itrlen] |
| PHCRID | `CRID` (ParticleID, algorithmType 3) | type = `ctlword`, PDG = best of (e, μ, π, K, p) by log-likelihood, likelihood = that value, parameters = [llik_e, llik_mu, llik_pi, llik_k, llik_p, liq_besthyp, gas_besthyp, liq_nhits, gas_nhits] |
| MCHEAD | event parameters `mc_ntot`, `mc_origin`, `mc_primary_quark_pdg`, `mc_ip_mm` | `origin` is a flag word: 1 uū, 2 dd̄, 4 ss̄, 8 cc̄, 16 bb̄ (only these five occur) |
| MCPART | `MCParticles` + `MCParticlesSLDType`, `MCParticlesSLDOrigin`, `MCParticlesSLDId` (int32) | see §2.3 |
| PHCHRG.imc | `ChargedTracksMCLinks` (Track ↔ MCParticle) | \|imc\| is the MCPART bank id; weight +1 for imc>0, −1 for imc<0 (sign meaning unknown, preserved) |
| MCPNT | `ReconstructedParticlesMCLinks` (RP ↔ MCParticle) | weight = `econtrib`; `reason`, `nhits` not converted |
| PHKCHRG | — | track–cluster matching quantities; the association itself is in PHPOINT |

## 2. Conventions

### 2.1 Tracks

`hlxpar = (φ, 1/pT, tanλ, x, y, z)` at the point (x, y, z) on the helix, which for
almost all tracks is the point of closest approach to the beam spot (|x,y| median
1.8 mm, same as the beam position). `1/pT` is unsigned; the charge is the separate
`charge` field. The `TrackState` is written at that point:

- `referencePoint = (x, y, z)·10`, `D0 = Z0 = 0`, `phi = φ`, `tanLambda = tanλ`,
- `omega = charge · 2.99792458·10⁻⁴ · 0.6 T / pT` mm⁻¹ (LCIO sign convention, B along +z),
- covariance: `dhlxpar` is the packed lower triangle of the 5×5 matrix of
  (φ, 1/pT, tanλ, TSI, ETA), where TSI is the displacement along the trajectory normal in
  the bending plane and ETA along the binormal (jazelle_reader). TSI is mapped to d0 and
  ETA to z0 (both ×10 for mm), 1/pT to omega (× charge·2.998e-4·B), and the terms are
  re-packed into EDM4hep's (d0, φ, ω, z0, tanλ, t) lower triangle. Median resolutions read
  back: σφ 1.7–2.3 mrad, σd0 50–130 µm, σz0 50–130 µm.

`tkpar0`, `tkpar[5]`, `dtkpar[15]` (a second parameter set, probably at the calorimeter
face) are not converted: their meaning is not documented in the reader.

### 2.2 Clusters

PHKLUS gives only angles. The position is placed on the nominal LAC front face: barrel
radius 1768 mm for |cosθ| ≤ 0.84, endcap plane |z| = 2740 mm beyond. These two numbers
are nominal; the exact angles are kept in `shapeParameters`, so the position can be
recomputed. `status` values are 2097152 (1≪21) plus low bits; the meaning of the bits is
not documented and the word is kept in `type`.

### 2.3 Monte Carlo particles

**Particle codes.** `ptype` uses the JETSET 6.3 numbering for |code| < 1000 (7 e⁻, 8 νₑ,
9 μ⁻, 10 ν_μ, 11 τ⁻, 12 ν_τ; 17 π⁺, 18 K⁺, 19 K⁰, 20–22 D, 23 π⁰, 24 η, 25 η′, 26 η_c,
27–29 ρ⁺/K*, 30–32 D*, 33 ρ⁰, 34 ω, 35 φ, 36 J/ψ, 37 K_S, 38 K_L; 41 p, 42 n, 43–47 Σ/Ξ,
48–51 Σ_c/Ξ_c′, 57 Λ, 58 Λ_c, 59–60 Ξ_c, 61–64 Δ, 65–69 Σ*/Ξ*, 70 Ω⁻, 71–75 Σ_c*/Ξ_c*;
partons 500 g, 501 u, 502 d, 503 s, 504 c, 505 b in the JETSET 6.x flavour order) with
SLD extensions (6/16 a₁, 94/97 K₁, 101–103 B, 123–125 B*, 145–147 Σ_b, 241–243 Λ_b/Ξ_b,
293–297 Σ_b*/Ξ_b*, 91–93 d/t/α from GEANT), and `1000000 + KF` in the JETSET 7.3
numbering above (1000091 cluster, 1000092 string, 1000094 CM-shower system, 1000445 χ_c2,
1010221 f₀(980), 1010321 K₀*(1430), 1010441 χ_c0, 1020443 χ_c1, 1030443 ψ(2S)).
**A positive code is the particle, a negative code the antiparticle.** The identifications
were made from the median mass of the stored energy, the charge, and the decay products
(e.g. Σ⁰ → 1 + 57, Λ → 41 + 17, Δ⁺⁺ → 41 + 17), and the whole convention is checked by
the net charge of the generator-level final state, which is exactly zero in every event.
The complete table with masses and notes is `sld_edm4hep/ptypes.py`.

**Charge and mass.** The `charge` column of the Parquet export is not used: the reader's
decoding has the sign inverted and the 1/3 ↔ 2/3 magnitudes swapped (proton +41 → −1,
deuteron 91 → −1, c quark 504 → −1/3). Charge comes from the table. The `e` column is
also not used for known hadrons and leptons: the bank omits the energy word for many
codes (γ, π, K, K⁰, K_S, K_L, p, n, ν, partons, …) and the reader then inserts a mass
from a table that is wrong for this corpus (7–10 are leptons, not K/p/n; 17 is π, not
6), so the Parquet energy of every pion is |p| and of every electron is √(p²+m_K²).
The EDM4hep mass is the tabulated one; only partons and colour-singlet systems
(500–505, 1000091/92/94) take √(E²−p²) from the stored energy. Twelve excited-charm
codes (109–120, a few per 10⁴ events) and six objects of mass 3.4–3.8 GeV decaying to
Λ_c + (e, γ, Λ) (1000073–78, 17 per 60k events) have no PDG identification; they carry
PDG 0 and their SLD code in `MCParticlesSLDType`.

**Tree.** `parent_id` is a bank id, resolved through the `id` column (ids are not
contiguous because the mini-DST keeps 243 of the ~885 generated entries per event, but
parents always precede daughters); 99.2% of entries have a parent, the rest (parent_id 0,
the primary quarks and ISR photons) are roots. `xt` is the termination/decay point and
becomes `endpoint`; `vertex` is the parent's endpoint, or the MCHEAD interaction point
for roots.

**Generator vs. simulation.** The bank contains the generator record and, for each
particle tracked through the detector, a GEANT *copy* whose parent is the generator entry
of the same particle (bit SWUM, 1≪20). A particle is generator-level if it is a root or if
its parent is generator-level and decayed in the generator (bit DECAYED, 1≪0, and not
DECAYFLT 1≪1 or INTERACT 1≪5). Generator-level particles get `generatorStatus` 1 (final),
2 (decayed), 21 (quarks, gluons, strings, clusters; documented choice), 4 (BEAM bit);
everything else gets 0 with `simulatorStatus` bit CreatedInSimulation. The STOPPED bit
(1≪4) sets the Stopped bit. On the released sample the CreatedInSimulation fraction equals
the SWUM fraction (62.0%) and a hadronic Z event has 41.7 generator-level final-state
particles. The raw origin word is kept in `MCParticlesSLDOrigin`; its other bits
(DECAYFLT, BEAMPIPE, NOINTER, STOPPED, INTERACT, INTSHDEP, PRIMARY, ISTOP1/2, KALTOCDC,
SWERROR, SW2MNYST, SWOUTOFT, EMAXTERM, NOTTRACK, ISR, BEAM, PREFRAG, SWUM) follow the
jazelle_reader `MCPART::Origin` enumeration.

## 3. Provenance

Each output file carries a `metadata` frame with the input file name and its SHA-256,
`sld_edm4hep_version`, `mapping_version`, units, the field orderings of
`subdetectorHitNumbers` and `shapeParameters`, the ParticleID algorithm ids and the
complete code → PDG table used.

## 4. Not resolved

1. `PHCHRG.dedx`: a packed 32-bit word (e.g. 318742740); the packing is not in the
   reader or the corpus. Kept raw in `ChargedTracksDedxRaw`; `ChargedTracksDqDx` is not
   produced rather than filled with a guess.
2. `PHCHRG.tkpar/dtkpar`, `smwstat`, `vxdhit`, `mustat`, `estat`, `chi2`, `chi2v`: not
   converted (meaning or reference surface not documented).
3. Sign of `PHCHRG.imc` (2.4% of tracks have imc ≤ 0): preserved as the link weight sign; imc = 0 gives no link.
4. LAC front-face radius and endcap z: nominal values; angles preserved.
5. Cluster and track `status` bit meanings: words kept raw in `type`.
6. Codes 109–120 and 1000073–78: no PDG assignment (see §2.3).
7. The muon sample asymmetry (SLD code −9 four times more frequent than +9 among GEANT
   secondaries) is a property of the input, not of the conversion; noted here so that
   nobody attributes it to the code/sign convention, which the net-charge test fixes.
