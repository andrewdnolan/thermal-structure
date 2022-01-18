################################################################################
# Makefile struggles with rules which generate multiple targets. We use pattern
# rule to make sure the rule is not run once for each target (i.e. generated file).
# See:
# 		https://stackoverflow.com/a/3077254/10221482
################################################################################

template = run/group.in run/group.out run/group_spline_init.sh

small  = $(subst group,small,$(template))
medium = $(subst group,medium,$(template))
large  = $(subst group,large,$(template))

all: $(small) $(medium) $(large)

$(subst .,%,$(small)): prepare2submit
	./prepare2submit --force --time="1:15:00" --mem="200MB" --group="small"

$(subst .,%,$(medium)): prepare2submit
	./prepare2submit --force --time="2:30:00" --mem="500MB" --group="medium"

$(subst .,%,$(large)): prepare2submit
	./prepare2submit --force --time="3:30:00" --mem="1GB" --group="large"
