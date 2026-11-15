#5/2/2024
import six
from cobra import Reaction, Model

def mod2irrev(model):
	rxnid=[]
	rxnname=[]
	rgpr=[]
	subsys=[]
	lbs=[];ubs=[];rmets=[];revs=[]
	for r in model.reactions:
			# If a reaction is reverse only, the forward reaction (which
			# will be constrained to 0) will be left in the model.
			if r.lower_bound < 0:
					rxnid.append(r.id+'_f')
					rxnname.append(r.name)
					rgpr.append(r.gpr)
					subsys.append(r.subsystem)
					lbs.append(0)
					ubs.append(max(0,r.upper_bound))
					rxn_dict = {k: v  for k, v in six.iteritems(r._metabolites)}
					rmets.append(rxn_dict)
					#backword
					rxnid.append(r.id+'_b')
					rxnname.append(r.name)
					rgpr.append(r.gpr)
					subsys.append(r.subsystem)
					lbs.append(0)
					ubs.append(r.lower_bound*-1)
					rxn_dict = {k: v * -1 for k, v in six.iteritems(r._metabolites)}
					rmets.append(rxn_dict)
			else:
					rxnid.append(r.id)
					rxnname.append(r.name)
					rgpr.append(r.gpr)
					subsys.append(r.subsystem)
					lbs.append(r.lower_bound)
					ubs.append(r.upper_bound)
					rxn_dict = {k: v  for k, v in six.iteritems(r._metabolites)}
					rmets.append(rxn_dict)
					revs.append(r.reversibility)
	
	nr=len(rxnid)
	print(nr)
	rl=[]
	for i in range(nr):
		rxn= Reaction(rxnid[i])
		rxn.name=rxnname[i]
		rxn.gpr=rgpr[i];rxn.subsytem=subsys[i]
		rxn.lower_bound=lbs[i];rxn.upper_bound=ubs[i]
		rxn.add_metabolites(rmets[i])
        rxn.reversibility=False
		rl.append(rxn)
	
	irm = Model('irrev:'+str(model.name))
	irm.add_reactions(rl)
	irm.objective=model.objective
    
	return(irm)