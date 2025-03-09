function C = Npc(p)
% C = NPC(p); generates comparison numbers of photo-counts per channel
% assumes output unchanged from input, apart from transmission factor t
% eg, identity matrix or uniform thermal+unitary

C = (sinh(p.sqz).^2+p.alpha.^2).*p.tr.^2;
C = reshape(C,[1,1,numel(C)]);
end                                         %End nc function