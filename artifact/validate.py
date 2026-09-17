"""Meaningful model checks; no assertion that a modality must win or lose."""
import json
import math
from pathlib import Path
import numpy as np
from cryogenic_model import (CONFIG_PATH, added_wall_energy, case_result,
    cooling_factor, electrical_multiplier, energy_gain, reference_score, shot_count)


def dense_reference(angles):
    """Independent tensor-axis simulator, with 2-qubit CNOT matrix contraction."""
    state=np.zeros([2]*8,dtype=complex); state[(0,)*8]=1
    def apply(gate,qubits):
        nonlocal state
        axes=[7-q for q in qubits]
        rest=[axis for axis in range(8) if axis not in axes]
        order=axes+rest
        moved=np.transpose(state,order).reshape(2**len(axes),-1)
        state=np.transpose((gate@moved).reshape([2]*8),np.argsort(order))
    def ry(a):
        return np.array([[np.cos(a/2),-np.sin(a/2)],[np.sin(a/2),np.cos(a/2)]])
    for q,a in enumerate(angles): apply(ry(a),[q])
    cx=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    for l in range(4):
        for q in range(8):
            apply(ry(np.sin(1+8*l+q)/3),[q])
            a=np.cos(1+8*l+q)/3
            apply(np.diag(np.exp(np.array([-1j,1j])*a/2)),[q])
        for q in range(8): apply(cx,[q,(q+1)%8])
    probs=np.abs(state)**2
    return float(probs[...,0].sum()-probs[...,1].sum())


def run():
    current=json.loads(CONFIG_PATH.read_text())
    case_result(current)  # Check the active scenario independently of golden fixtures.
    c=dict(current)
    c.update(hot_temperature_k=300., controller_temperature_k=4., carnot_fraction=.02,
             qubits=8, controller_power_per_qubit_w=.023, shot_duration_s=.0001,
             failure_probability=.05, statistical_error=.1, tokens_per_decision=32,
             availability=.9, token_rate_per_s=1000., available_controller_cooling_w_per_plant=1.5,
             fixed_wall_power_w_per_plant=5000.)
    assert cooling_factor(300)==0 and electrical_multiplier(300)==1
    assert electrical_multiplier(299.9)>1
    for t,w in [(4,3700),(1,14950),(.02,749950)]:
        assert math.isclose(cooling_factor(t),w)
    # A locally powered cold load requires direct input plus cooling.
    assert math.isclose(added_wall_energy(.01,[(4,.02,.01)]),37.01)
    # Warm input partly absorbed cold: do not add the absorbed energy twice.
    assert math.isclose(added_wall_energy(.1,[(4,.02,.001)]),3.8)
    assert shot_count(.2)==185 and shot_count(.1)==738 and shot_count(.05)==2952
    for eps in (.2,.1,.05):
        s=shot_count(eps)
        assert 2*math.exp(-s*eps**2/2)<=.05
        assert 2*math.exp(-(s-1)*eps**2/2)>.05
    f,h,a=.3,.05,3701
    d=a/(1-h/f)
    assert math.isclose(float(energy_gain(d,a,f,h)),1.0)
    assert energy_gain(d*.99,a,f,h)<1<energy_gain(d*1.01,a,f,h)
    r=case_result(c)
    assert math.isclose(r["controller_energy_j_per_decision"],50.2566192)
    assert r["engines_required"]==3 and r["engines_per_plant"]==8
    assert r["engines_required"]*c["availability"] >= c["token_rate_per_s"]*r["decision_duration_s"]/c["tokens_per_decision"]
    assert (r["engines_required"]-1)*c["availability"] < c["token_rate_per_s"]*r["decision_duration_s"]/c["tokens_per_decision"]
    assert case_result(c,rate=10000)["plants_required"]==4
    assert math.isclose(case_result(c,rate=10)["controller_energy_j_per_token"],r["controller_energy_j_per_token"])
    s=reference_score()
    assert abs(s["norm"]-1)<1e-12 and -1<=s["ideal_z0"]<=1
    # Independent dense tensor contraction checks the index-loop reference simulator.
    for offset in (-.2, 0, .2):
        angles=np.linspace(-.6,.8,8)+offset
        z=dense_reference(angles)
        assert math.isclose(reference_score(angles)["ideal_z0"],z,abs_tol=1e-12)
    # Capacity transitions: fixed allocations recur rather than vanishing at scale.
    limit=8*c["availability"]*c["tokens_per_decision"]/r["decision_duration_s"]
    assert case_result(c,rate=limit*(1-1e-10))["plants_required"]==1
    assert case_result(c,rate=limit*(1+1e-10))["plants_required"]==2
    assert case_result(c,eta=.1)["controller_energy_j_per_decision"] < r["controller_energy_j_per_decision"]
    try:
        reference_score([0]*7)
    except ValueError:
        pass
    else:
        raise AssertionError("reference accepted seven angles")
    # Invalid physical domains must be rejected.
    for args in [(0,.02),(301,.02),(4,0),(4,1.01),(float("nan"),.02)]:
        try:
            cooling_factor(*args)
        except ValueError:
            pass
        else:
            raise AssertionError(args)
    print("PASS: energy accounting, shot guarantee, crossover, capacity, reference circuit and input domains")


if __name__=="__main__":
    run()
