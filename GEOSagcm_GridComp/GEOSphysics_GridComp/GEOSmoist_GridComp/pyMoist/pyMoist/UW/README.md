
- UW Shallow Water Convection physics port - from the Moist component.

- This ports is aimed at the core numerics, e.g.: uwshcu.F90

- Version ported: v2.5.2 (as part of GEOS v11..5.2)

- Some errors still remain when testing UW in debug backend. I believe they are stemming from the
  'qtu' variable in the iter_xc loop in the buoyancy_sorting stencil. The values start off small,
  around 1e-10. Then they grow and spread, resulting in translate test failure. For now, an override
  has been put into place in the multimodal translate test so that UW passes. But we will want to
  revist this and solve these errors at some point.

- Last updated: 7/17/25
