# Reproductions — script → paper map

Auto-derived from `papers/*.tex` citations. `paperNN_*/` modules hold each paper's
reproduction drivers (run with `nwt-repro run <name>`); `supporting/` holds shared /
exploratory analysis scripts not tied to one paper's reproduction.


**85 paper-reproduction drivers across 11 modules; 81 supporting scripts.**


## Paper reproductions

| module | script | cited by |
|---|---|---|
| `paper06_form_factor` | `paper6_figures.py` | Paper 6 |
| `paper07_knot_invariants` | `nwt_charge_from_framing.py` | Paper 7 |
| `paper07_knot_invariants` | `nwt_jones_polynomial.py` | Paper 7 |
| `paper07_knot_invariants` | `nwt_surgery_network.py` | Paper 7 |
| `paper08_ewk_couplings` | `nwt_crossing_geometry.py` | Paper 8, Paper 9 |
| `paper08_ewk_couplings` | `nwt_ewk_boson_scan.py` | Paper 8, Paper 9 |
| `paper08_ewk_couplings` | `nwt_hopf_component_analysis.py` | Paper 8 |
| `paper08_ewk_couplings` | `nwt_knot_eigensolver_v3.py` | Paper 8, Paper 9 |
| `paper08_ewk_couplings` | `nwt_reidemeister_couplings.py` | Paper 8, Paper 9 |
| `paper09_lepton_observables` | `nwt_g_minus_2.py` | Paper 9 |
| `paper09_lepton_observables` | `nwt_gpe_knots_3d.py` | Paper 9 |
| `paper09_lepton_observables` | `nwt_pmns_3d.py` | Paper 9 |
| `paper13_spectrum` | `nwt_sensitivity_analysis.py` | Paper 13, Paper 21 |
| `paper17_qec` | `nwt_qec_KN_generalization.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_bit_quantum_from_bremermann.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_bracket_test.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_entanglement_structure.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_forward_prediction.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_heron_zne.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_interpretation_b_test.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_proportionality_constant.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_psl27_edge_transitivity.py` | Paper 17 |
| `paper17_qec` | `nwt_qec_route_a_so7_lift.py` | Paper 17 |
| `paper17_qec` | `nwt_spin7_fusion_compute.py` | Paper 17 |
| `paper17_qec` | `nwt_truncation_channel_mixing.py` | Paper 17 |
| `paper17_qec` | `nwt_truncation_mechanism.py` | Paper 17 |
| `paper17_qec` | `nwt_truncation_qdef.py` | Paper 17 |
| `paper18_gravity` | `nwt_paper18_G1_setup.py` | Paper 18 |
| `paper18_gravity` | `nwt_paper18_G2_matter_loop.py` | Paper 18 |
| `paper18_gravity` | `nwt_paper18_G3_K7_sakharov.py` | Paper 18 |
| `paper18_gravity` | `nwt_paper18_G4_linearized_einstein.py` | Paper 18 |
| `paper18_gravity` | `nwt_paper18_G5_curved_sakharov.py` | Paper 18 |
| `paper18_gravity` | `nwt_paper18_G6_einstein_variation.py` | Paper 18 |
| `paper19_substrate_figures` | `nwt_multiview_demo.py` | Paper 19 |
| `paper19_substrate_figures` | `nwt_muon_decay_heron_demo.py` | Paper 19 |
| `paper19_substrate_figures` | `nwt_paper19_V1_vacuum_energy.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_V2_condensate_scale.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_V4_bit_creation.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W12_holographic_K7.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_2_fcc_phonons.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_2b_k7_internal.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_2c2_floppy_structure.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_2c3_gauge_field.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_2c_central_force.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_2d_two_parameter.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3_gluing_scan.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3b_edge_gluing.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3c_psl27.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3d_S4_decomp.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3e_eg_matter_coupling.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3f_constructed_eg.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3g_D4_subdecomp.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3h_paper6_consolidation.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3i_neutrino_phase_soliton.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3j_neutrino_NLO_NNLO.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_3k_PMNS_from_K8.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_4_lorentz_alpha.py` | — |
| `paper19_substrate_figures` | `nwt_paper19_W3_phonon_dispersion.py` | — |
| `paper19_substrate_figures` | `paper18_torus_unification.py` | Paper 19 |
| `paper19_substrate_figures` | `paper19_fig_qed_compton.py` | Paper 19 |
| `paper19_substrate_figures` | `paper19_substrate_figures.py` | Paper 19 |
| `paper20_neutrinos` | `g2_z3_subgroups.py` | Paper 20 |
| `paper20_neutrinos` | `k7_L3_full_second_order.py` | Paper 20 |
| `paper20_neutrinos` | `k7_charged_lepton_U_ell.py` | Paper 20 |
| `paper20_neutrinos` | `k7_delta_CP_topological.py` | Paper 20 |
| `paper20_neutrinos` | `k7_orbit_wilson_seesaw.py` | Paper 20 |
| `paper20_neutrinos` | `k7_pmns_mass_hierarchy_test.py` | Paper 20 |
| `paper20_neutrinos` | `k7_pmns_nlo_corrections.py` | Paper 20 |
| `paper20_neutrinos` | `k7_pmns_nlo_derivation.py` | — |
| `paper20_neutrinos` | `k7_seesaw_gauge_invariance.py` | Paper 20 |
| `paper20_neutrinos` | `k7_sterile_neutrino_phenomenology.py` | Paper 20 |
| `paper20_neutrinos` | `k7_triality_absolute_sterile_masses.py` | Paper 20 |
| `paper20_neutrinos` | `k7_triality_delta_CP.py` | Paper 20 |
| `paper20_neutrinos` | `k7_triality_seesaw_verification.py` | Paper 20 |
| `paper20_neutrinos` | `k8_extension_orbit_decomposition.py` | Paper 20 |
| `paper21_boson_rep` | `exp12_phase1_simulator_pilot.py` | — |
| `paper21_boson_rep` | `exp12_phase2_destructive.py` | — |
| `paper21_boson_rep` | `exp12_phase2_hardware.py` | Paper 21 |
| `paper21_boson_rep` | `exp12_scrambled_mapping_control.py` | Paper 21 |
| `paper21_boson_rep` | `exp15_cross_sector_equivalence.py` | Paper 21 |
| `paper21_boson_rep` | `exp28_phase28e_braket.py` | — |
| `paper21_boson_rep` | `k7_sector_rule_lookelsewhere.py` | Paper 21 |
| `paper21_boson_rep` | `nwt_assignment_degeneracy.py` | Paper 21 |
| `paper21_boson_rep` | `nwt_complete_spectrum_paper13.py` | Paper 21 |
| `paper22_baryogenesis` | `nwt_baryogenesis_chirality.py` | Paper 22 |

## Supporting scripts (`supporting/`)

Shared/exploratory analysis not a single paper's reproduction:

- `nwt_2T_8dim_explicit_b2_8.py`
- `nwt_2T_8dim_lambda2_b2_6.py`
- `nwt_2T_character_7dim_b2_4.py`
- `nwt_2T_fusion_21channels_b2_4.py`
- `nwt_2T_in_2I_check_b2_3a.py`
- `nwt_2T_spin7_clifford_b2_12.py`
- `nwt_2T_spin7_full_b2_10.py`
- `nwt_2T_spin7_normalizer_b2_11.py`
- `nwt_2T_spin7_search_b2_9.py`
- `nwt_87_prefactor_b2_14.py`
- `nwt_alpha_ladder.py`
- `nwt_beta_decay_landauer.py`
- `nwt_compressibility_scalar_b2_2.py`
- `nwt_conservation_laws.py`
- `nwt_crossing_phase.py`
- `nwt_emergent_c.py`
- `nwt_emergent_yang_mills.py`
- `nwt_eulerian_amplitude_b2_7.py`
- `nwt_ewk_decay_check.py`
- `nwt_fano_crossings_b2_1c.py`
- `nwt_gravity_structure_b2_1.py`
- `nwt_hyperon_landauer_survey.py`
- `nwt_k7_so7_wilson_b2_13.py`
- `nwt_knot_eigensolver.py`
- `nwt_knot_eigensolver_v2.py`
- `nwt_lagrangian_L1_fields.py`
- `nwt_lagrangian_L1b_uv_completion.py`
- `nwt_lagrangian_L2_kinetic_bps.py`
- `nwt_lagrangian_L3_skyrme_hopf.py`
- `nwt_lagrangian_L4_paper6_mass_spectrum.py`
- `nwt_lagrangian_L5_gravity_hierarchy.py`
- `nwt_lifetime_geometry.py`
- `nwt_mass_from_abelian_higgs.py`
- `nwt_multimode_carrier_scan.py`
- `nwt_neutrino_2sqrtalpha_origin.py`
- `nwt_neutrino_alpha6_decomposition.py`
- `nwt_neutrino_hierarchy_ladder.py`
- `nwt_neutrino_mass_seesaw.py`
- `nwt_nlo_alpha7_b2_15.py`
- `nwt_null_worldtube_kappa.py`
- `nwt_pmns_cross_topology.py`
- `nwt_pmns_nonlinear.py`
- `nwt_poincare_sphere_b2_0.py`
- `nwt_qec_bps_compton_bridge.py`
- `nwt_qec_heron_KN.py`
- `nwt_qec_heron_exp4.py`
- `nwt_qec_heron_exp5.py`
- `nwt_qec_heron_experiment.py`
- `nwt_qec_heron_fetch.py`
- `nwt_qec_syndrome_attractor.py`
- `nwt_qec_time_evolution.py`
- `nwt_qec_zne_continuation.py`
- `nwt_qec_zne_ratio.py`
- `nwt_qec_zne_reanalysis.py`
- `nwt_rational_decomposition_methodology.py`
- `nwt_residual_histogram.py`
- `nwt_rybakov_path_a.py`
- `nwt_self_consistent_knot_v3.py`
- `nwt_spin7_chain_b2_5.py`
- `nwt_surface_balance.py`
- `nwt_trefoil_orbit_linking_b2_3.py`
- `nwt_tt_scalar_resolution_b2_1a.py`
- `nwt_universality_argument.py`
- `nwt_volovik_bogoliubov.py`
- `nwt_volovik_c.py`
- `nwt_volovik_closure.py`
- `nwt_volovik_part_b.py`
- `nwt_volovik_two_mode.py`
- `nwt_vortex_fluctuations_b1_0.py`
- `nwt_vortex_fluctuations_b1_1.py`
- `nwt_vortex_fluctuations_b1_2.py`
- `nwt_vortex_fluctuations_b1_3.py`
- `nwt_vortex_fluctuations_b1_4.py`
- `nwt_vortex_fluctuations_b1_5.py`
- `nwt_vortex_gravity_flat.py`
- `nwt_zeta_phase0_scaffold.py`
- `nwt_zeta_phase1_free_scalar.py`
- `nwt_zeta_phase2_trefoil_bps.py`
- `nwt_zeta_phase3_trefoil_casimir.py`
- `nwt_zeta_phase4_curvature_corrections.py`
- `nwt_zeta_phase5_1overG.py`
