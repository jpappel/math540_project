all: proposal/proposal.pdf

proposal/proposal.pdf: proposal/proposal.md
	pandoc -o $@ $<
