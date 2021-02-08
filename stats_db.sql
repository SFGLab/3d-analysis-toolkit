-- HG00512
select sum(case when HG00512_anchor = 1 and HG00512_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when HG00512_anchor = 0 and HG00512_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- HG00513
select sum(case when HG00513_anchor = 1 and HG00513_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when HG00513_anchor = 0 and HG00513_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- HG00514
select sum(case when HG00514_anchor = 1 and HG00514_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when HG00514_anchor = 0 and HG00514_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- HG00731
select sum(case when HG00731_anchor = 1 and HG00731_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when HG00731_anchor = 0 and HG00731_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- HG00732
select sum(case when HG00732_anchor = 1 and HG00732_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when HG00732_anchor = 0 and HG00732_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- HG00733
select sum(case when HG00733_anchor = 1 and HG00733_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when HG00733_anchor = 0 and HG00733_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- NA19238
select sum(case when NA19238_anchor = 1 and NA19238_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when NA19238_anchor = 0 and NA19238_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- NA19239
select sum(case when NA19239_anchor = 1 and NA19239_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when NA19239_anchor = 0 and NA19239_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 
union all
-- NA19240
select sum(case when NA19240_anchor = 1 and NA19240_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "LOOP",
sum(case when NA19240_anchor = 0 and NA19240_sv = "0/1" and "SV_type" = "<INV>" then 1 else 0 end) as "NO LOOP"
from anchor_database 

