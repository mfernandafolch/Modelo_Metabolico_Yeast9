import warnings
import cobra


def _first_existing_met_ids(model, candidates):
    """Return the subset of candidate metabolite IDs that exist in the model."""
    existing = []
    for mid in candidates:
        try:
            model.metabolites.get_by_id(mid)
            existing.append(mid)
        except KeyError:
            pass
    return existing


def _remove_mets_from_reaction(model, rxn_id, met_ids):
    """
    Remove given metabolites from a reaction.
    """
    try:
        rxn = model.reactions.get_by_id(rxn_id)
    except KeyError:
        warnings.warn(f"Reaction '{rxn_id}' not found.")
        return [], met_ids

    removed = []
    missing = []

    for mid in met_ids:
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            missing.append(mid)
            continue

        coeff = rxn.metabolites.get(met, 0.0)
        if coeff != 0.0:
            rxn.add_metabolites({met: -coeff})
            removed.append(mid)

    return removed, missing


def modify_yeast9_anaerobic(model, in_place=False, verbose=True):
    """
    Applies the MATLAB modifications to Yeast9 model.
    """

    if not in_place:
        model = model.copy()

    # --- Remove metabolites from biomass reaction r_4598 ---
    mets_with_comp = [
        "s_3714[c]",
        "s_1198[c]",
        "s_1203[c]",
        "s_1207[c]",
        "s_1212[c]",
        "s_0529[c]",
    ]

    existing = _first_existing_met_ids(model, mets_with_comp)

    if len(existing) == len(mets_with_comp):
        mets_to_remove = mets_with_comp
    else:
        mets_to_remove = [
            "s_3714",
            "s_1198",
            "s_1203",
            "s_1207",
            "s_1212",
            "s_0529",
        ]

    removed, missing = _remove_mets_from_reaction(
        model, "r_4598", mets_to_remove
    )

    if verbose:
        print("Removed from r_4598:", removed)
        print("Not found:", missing)

    # --- Block pathways ---

    def set_lb(rid, value):
        try:
            model.reactions.get_by_id(rid).lower_bound = value
        except KeyError:
            if verbose:
                print(f"Reaction {rid} not found.")

    def set_ub(rid, value):
        try:
            model.reactions.get_by_id(rid).upper_bound = value
        except KeyError:
            if verbose:
                print(f"Reaction {rid} not found.")

    # Block oxaloacetate-malate shuttle
    set_lb("r_0713", 0.0)
    set_lb("r_0714", 0.0)

    # Block glycerol dehydrogenase
    set_ub("r_0487", 0.0)

    return model