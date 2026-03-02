# modificacion_yeast9v2.py

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
    """Remove given metabolites from a reaction (set their coefficient to 0)."""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
    except KeyError:
        warnings.warn(f"Reaction '{rxn_id}' not found.")
        return [], list(met_ids)

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
            rxn.add_metabolites({met: -coeff})  # cancels it to 0
            removed.append(mid)

    return removed, missing


def modify_yeast9_change2_only(model, in_place=False, verbose=True):
    """
    Applies ONLY the MATLAB 2nd change to Yeast9:
    Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
    by setting their stoichiometric coefficients to 0 in reaction r_4598.

    Returns the modified model (copy by default).
    """
    if not in_place:
        model = model.copy()

    # Try with compartment tags first (like MATLAB)
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
        # Fallback without compartment tags
        mets_to_remove = [
            "s_3714",
            "s_1198",
            "s_1203",
            "s_1207",
            "s_1212",
            "s_0529",
        ]

    removed, missing = _remove_mets_from_reaction(model, "r_4598", mets_to_remove)

    if verbose:
        print("Removed from r_4598:", removed)
        print("Not found:", missing)

    return model