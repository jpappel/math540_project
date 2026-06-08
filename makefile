
# presentation related variables
PIC_SRC := $(wildcard presentation/pictures/*.jl)
PICS := $(PIC_SRC:.jl=.png)
ANIM_SRC := $(wildcard presentation/animations/*.jl)
ANIMS := $(ANIM_SRC:.jl=.mp4)

.PHONY: all info

all: proposal/proposal.pdf presentation/presentation.html

proposal/proposal.pdf: proposal/proposal.md
	pandoc -o $@ $<

presentation/presentation.html: presentation/presentation.md $(PICS) $(ANIMS)
	pandoc -t revealjs \
		-s \
		--slide-level=3 \
		-V revealjs-url="https://unpkg.com/reveal.js@^4" \
		-o $@ $<

presentation/pictures/%.png: presentation/pictures/%.jl
	julia $<
	mv $(notdir $@) $@

presentation/animations/%.mp4: presentation/animations/%.jl
	julia $<
	mv $(notdir $@) $@

info:
	@echo PIC_SRC: $(PIC_SRC)
	@echo PICS: $(PICS)
	@echo ANIM_SRC: $(ANIM_SRC)
	@echo ANIMS: $(ANIMS)
