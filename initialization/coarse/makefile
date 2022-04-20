################################################################################
# Makefile struggles with rules which generate multiple targets. We use pattern
# rule to make sure the rule is not run once for each target (i.e. generated file).
# See:
# 		https://stackoverflow.com/a/3077254/10221482
################################################################################

out_temp = run/group.in run/group.out run/group_spline_init.sh
param_temp = params/glac.json

small_out = $(subst group,small,$(out_temp))
small_glacs = crmpt12 crmpt18-a crmpt18-b glc1-a glc1-b
small_params = $(foreach glac,$(small_glacs),$(subst glac,$(glac),$(param_temp)))

medium_out= $(subst group,medium,$(out_temp))
medium_glacs = lilk-a lilk-b klun-a klun-b sprg
medium_params = $(foreach glac,$(medium_glacs),$(subst glac,$(glac),$(param_temp)))

large_out = $(subst group,large,$(out_temp))
large_glacs = twds-a twds-b fish klut-a klut-b
large_params = $(foreach glac,$(large_glacs),$(subst glac,$(glac),$(param_temp)))

all: $(small_out) $(medium_out) $(large_out)

$(subst .,%,$(small_out)): prepare2submit $(small_params)
	./prepare2submit --force --time="5:45:00" --mem="1250M" --group="small" -s 30

$(subst .,%,$(medium_out)): prepare2submit $(medium_params)
	./prepare2submit --force --time="9:30:00" --mem="2250M" --group="medium" -s 25

$(subst .,%,$(large_out)): prepare2submit $(large_params)
	./prepare2submit --force --time="10:00:00" --mem="2250M" --group="large" -s 25
